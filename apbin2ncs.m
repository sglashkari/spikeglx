%% apbin2ncs reads a Spike GLX AP bin files to Neuralynx CSC files
% It also remvoes the DC offset, does common average referencing (CAR) and
% bandpass filtering:
%   600-6000 Hz: CSC for AP
%   1-400 Hz: LFP (only for the last channel)
% At the end, it calls Vyash's code to extract the params and waveform
% files.
%
% *.ap.bin ==> *.ncs
%
%   See also APBIN2NCS, WINCLUST2MAT.
%
% Date 2022-04-29
% Author Shahin G Lashkari
%
% Check:
% https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/UserManual.md#output-file-format-and-tools
%
clc; clear
close all
numCh=32;
answer = inputdlg('Enter a value for threshold (in microVolts)','Threshold', [1 45],{'80'});
thresh = str2double(answer);
fprintf('Threshold is %d uV.\n',thresh)
%% Read data from files
[binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select One or More Binary Files','MultiSelect','on');
if isa(binaryFile,'double')
    return;
elseif isa(binaryFile,'char')
    binaryFile = {binaryFile};
end

[csvFile,csvPath] = uigetfile('D:\ChannelLog\*.csv','Channels of Interest: Select One or More CSV Files','MultiSelect','on');
if isa(csvFile,'double')
    return;
elseif isa(csvFile,'char')
    csvFile = {csvFile};
end

selparentpath = uigetdir('D:\Analysis','Select the main Directory for Saving CSC Files');
if isa(selparentpath,'double')
    return;
end

chunksize = 300; % 300 samples of length 512, almost 5 seconds: 5 sec x 30Khz / 512 sample

% updated 4/29 https://billkarsh.github.io/SpikeGLX/Sgl_help/Metadata_30.html
meta = ReadMeta(binaryFile{1}, path);
nChan = str2double(meta.nSavedChans); % 385
if meta.imDatPrb_type == "0"
    disp('Neuropixels 1.0 is selected!')
    NP = 1;
    peak2peak = 5*1e-3; % 5 mV for NP 1.0
    bits = 10; % 10-bit for NP 1.0
    % For type 0 imec probes:
    Imax = str2double(meta.imMaxInt);
    Vmax = str2double(meta.imAiRangeMax);
    % reading AP gains from IMRO table in the meta file
    C = textscan(meta.imroTbl, '%s', 'Delimiter', ')(' );
    C = {C{1}{4:2:end}};
    gain = zeros(1,nChan);
    for i=1:nChan-1
        imroArray = cell2mat(textscan(C{i}, '%d %d %d %d %d', 'Delimiter', ' '));
        gain(i) = imroArray(4);
    end
elseif meta.imDatPrb_type == "21" || meta.imDatPrb_type == "24"
    disp('Neuropixels 2.0 is selected!')
    NP = 2;
    peak2peak = 12.5*1e-3; % 12.5 mV for NP 2.0
    bits = 14; % 14-bit for NP 2.0
    % For type 21 or type 24 imec probes:
    Imax = str2double(meta.imMaxInt);
    Vmax = str2double(meta.imAiRangeMax);
    gain = 80*ones(1,nChan);
else
    error(['Type ' meta.imDatPrb_type ' is not supported!'])
end
% voltperbit = peak2peak/2^bits
voltperbit = Vmax ./ Imax ./ gain;
Nlx_ADBitVolts = 0.000000036621093749999997;
Nlx_bits_per_NP_bits = voltperbit/Nlx_ADBitVolts;


Vyashpath = 'C:\Users\neuropixels\Neuropixels\NPtoWinclust_Pipeline\makeParms';
Vyashpath = uigetdir(Vyashpath,'Select Vyash''s Code Directory');
if isa(Vyashpath,'double')
    return;
end
thisFilePath = pwd;

%%
start = tic;
disp('Starting ...')
for l = 1:length(csvFile)
    fprintf('\nChannel list %d from %d channel lists.\n\n',l,length(csvFile));
    
    range = readmatrix(fullfile(csvPath, csvFile{l})) + 1; % selection range in matlab is 1..385, while for the range for AP is 0..384
    range = reshape(range,1,[]); % making range a row vector
    selpath = fullfile(selparentpath, extractBefore(csvFile{l},'.csv'));
    
    if ~exist(selpath, 'dir')
        fprintf('Creating %s directory ...\n', selpath)
        mkdir(selpath);
        
        %% Read AP binary files and write them into csc files
        cd(thisFilePath);
        disp('Reading and writing started ...')
        start1 = tic;
        for k = 1:length(binaryFile)
            meta = ReadMeta(binaryFile{k}, path);
            nChan = str2double(meta.nSavedChans); % 385
            nSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
            n = floor(nSamp/512/chunksize);
            dt = 1/str2double(meta.imSampRate); % in seconds
            tStart{k} = posixtime(datetime(strrep(meta.fileCreateTime,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone', 'America/New_York'));
            
            %time = 0:dt:meta.fileTimeSecs;
            
            % opening
            binFile = fopen(fullfile(path, binaryFile{k}), 'rb');
            
            if k == 1
                cscFile = cell(nChan,1); % 385 x 1
                for j=range
                    FileName = [selpath filesep 'CSC' num2str(j-1) '.ncs'];
                    cscFile{j} = fopen(FileName,'wb');
                end
                
                % Header
                headerSize = 16*1024; % 16 kb
                dummy = fopen('header.txt','rb');
                hdr = fread(dummy, [1 headerSize], '*char');
                fclose(dummy);
                for j=range
                    fwrite(cscFile{j}, hdr, 'char'); % Write Header of CSCs
                end
                disp('Headers written to the files ...')
            end
            
            % Body
            for i=1:n-1
                Samples = fread(binFile, [nChan, 512*chunksize], 'int16=>double')'; % (512xchunksize) x 385
                Samples(:,range) = Samples(:,range) - mean(Samples(:,range),1); % remove DC offset
                Samples(:,range) = Samples(:,range) - median(Samples(:,range),2); % CAR (changed from mean to median SGL 2021-01-21)
                Samples(:,range) = - Samples(:,range); % INVERTED
                
                TimeStamp = (tStart{k} - tStart{1} + (i-1)*512*chunksize*dt:512*dt:(i*512*chunksize-1)*dt )*1e6; % in microseconds
                Ttime = (i-1)*512*chunksize*dt:dt:(i*512*chunksize-1)*dt;
               
                SampleFreq=str2double(meta.imSampRate);
                NumValidSamples=512;
                
                % Write Body of CSCs
                for j=range
                    ChannelNumber=j-1;
                    for l = 1:chunksize
                        fwrite(cscFile{j}, TimeStamp(l), 'uint64');
                        fwrite(cscFile{j}, ChannelNumber, 'uint32');
                        fwrite(cscFile{j}, SampleFreq, 'uint32');
                        fwrite(cscFile{j}, NumValidSamples, 'uint32');
                        % 512 x N
                        fwrite(cscFile{j}, Nlx_bits_per_NP_bits(j) * Samples((l-1)*512+1:l*512,j), 'int16'); % Convert NP bits to NLX bits
                    end
                end
                if mod(i,5*round(n/100))==0 % display progress every 5 percent
                    fprintf('Read/Write: %.0f seconds, %.0f%% done.\n',toc(start1), i/(n-1)*100)
                end
            end
            fclose(binFile);
        end
        % closing
        for j=range
            fclose(cscFile{j});
        end
        fprintf(['It took ' datestr(seconds(toc(start1)),'HH:MM:SS') ,' to read and write files.\n\n']);
        
        %% Bandpass filter: AP 600 - 6000 Hz, LFP 1 - 400 Hz
        disp('Filtering started ...')
        start2= tic;
        for j=range
            [time,data,header,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc([selpath filesep 'CSC' num2str(j-1) '.ncs']);
            data_filtered = filterlfp(time, data, 600, 6000);
            write_bin_csc([selpath filesep 'CSC' num2str(j-1) '.ncs'], time,data_filtered,header,ChannelNumber,SampleFreq,NumValidSamples);
            %if j==range(end)
                data_filtered = filterlfp(time, data, 1, 400);
                write_bin_csc([selpath filesep 'LFP' num2str(j-1) '.ncs'], time,data_filtered,header,ChannelNumber,SampleFreq,NumValidSamples);
            %end
            fprintf('\rFiltering: %.0f seconds, %.0f%% done.',toc(start2), find(range==j)/length(range)*100)
        end
        fprintf(['\nIt took ' datestr(seconds(toc(start2)),'HH:MM:SS') ,' to filter ', num2str(length(range)), ' csc files.\n\n']);
        
        %% Prep for Vyash's Analysis
        selpathB = [selpath filesep 'B'];
        mkdir(selpathB);
        ncsFiles = dir([selpath filesep '*.ncs']);
        i = 0;
        for j=range
            i = i+1;
            movefile([selpath filesep 'CSC' num2str(j-1) '.ncs'], [selpathB filesep 'CSC_B' num2str(i) '.ncs']);
            %if j==range(end)
                movefile([selpath filesep 'LFP' num2str(j-1) '.ncs'], [selpathB filesep 'LFP' num2str(i) '.ncs']);
            %end
        end
        copyfile('VideoReport', [selpathB filesep 'VideoReport']);
        copyfile('VideoSampling', [selpathB filesep 'VideoSampling']);
        
    else
        fprintf('%s directory already exists ...\n', selpath)
        selpathB = [selpath filesep 'B'];
    end
    %% Vyash's Analysis
    cd(Vyashpath);
    pause(1);
    start2= tic;
    
    if ~exist([selpathB filesep 'B_' num2str(thresh)], 'dir')
        system(['python makePolytrodeSpikes_singleChannel.py ' selpathB ' B ' num2str(numCh) ' ' num2str(thresh)]);
        system(['python expandCSCFile.py ' selpathB filesep 'B_' num2str(thresh)]);
        system(['python makeParms.py ' selpathB filesep 'B_' num2str(thresh)]);
        
        delete([selpathB filesep 'B_' num2str(thresh) filesep '*.npy']);
        delete([selpathB filesep 'B_' num2str(thresh) filesep '*.npz']);
        delete([selpathB filesep 'B_' num2str(thresh) filesep '*.memmap']);
        
        fprintf(['\nIt took ' datestr(seconds(toc(start2)),'HH:MM:SS') ,' to create the parms file.\n\n']);
    else
        fprintf('%s directory already exists ...\n', [selpathB filesep 'B_' num2str(thresh)])
    end
    beep;
end
%% Done
cd(thisFilePath);
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);