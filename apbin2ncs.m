%% apbin2ncs reads a Spike GLX AP bin files to Neuralynx CSC files
% It also remvoes the DC offset, does common average referencing (CAR) and
% bandpass filtering (600-6000 Hz).
%
% *.ap.bin ==> *.ncs
% This program is written by Shahin G Lashkari
% Last upgrade: 2021-11-13
%
% Check:
% https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/UserManual.md#output-file-format-and-tools
%
clc; clear
close all

%% Read data from files
[binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select One or More Binary Files','MultiSelect','on');
if isa(binaryFile,'double')
    return;
elseif isa(binaryFile,'char')
    binaryFile = {binaryFile};
end

[csvFile,csvPath] = uigetfile('*.csv','Channels of Interest: Select a CSV File to Open');
if isa(csvFile,'double')
    answer=questdlg('No channel is selected! Do you want to include all channels?', 'Warning!', 'Yes', 'No', 'No');
    if ~isequal(answer,'Yes')
        disp('No channel is selected!')
        return
    else
        range = 1:32;
        disp('First 32 channels are selected.')
    end
else
    csvFile = fullfile(csvPath, csvFile);
    range = readmatrix(csvFile) + 1; % selection range in matlab is 1..385, while for the range for AP is 0..384
end

selpath = uigetdir('D:\Analysis','Select a Directory for Saving CSC Files');
if isa(selpath,'double')
    return;
end

%% Read AP binary files and write them into csc files
answer=questdlg('Do you use Neuropixels 1.0 or 2.0?', 'Warning!', 'NP 1.0', 'NP 2.0', 'NP 2.0');
if isequal(answer,'NP 1.0')
    disp('Neuropixels 1.0 is selected!')
    NP = 1;
    peak2peak = 5*1e-3; % 5 mV for NP 1.0
    bits = 10; % 10-bit for NP 1.0
else
    disp('Neuropixels 2.0 is selected!')
    NP = 2;
    peak2peak = 12.5*1e-3; % 12.5 mV for NP 2.0
    bits = 14; % 14-bit for NP 2.0
end
chunksize = 300; % almost 5 seconds: 5 sec x 30Khz / 512 sample
Nlx_ADBitVolts = 0.000000036621093749999997;
voltperbit = peak2peak/2^bits/Nlx_ADBitVolts;

disp('Reading and writing started ...')
start = tic;
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
        Samples(:,range) = Samples(:,range) - mean(Samples(:,range),2); % CAR
        Samples = -voltperbit * Samples; % microvolts to bits INVERTED
        
        TimeStamp = (tStart{k} - tStart{1} + (i-1)*512*chunksize*dt:512*dt:(i*512*chunksize-1)*dt )*1e6; % in microseconds
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
                fwrite(cscFile{j}, Samples((l-1)*512+1:l*512,j), 'int16');
            end
        end
        if mod(i,5*round(n/100))==0 % display progress every 5 percent
            fprintf('Read/Write: %.0f seconds, %.0f%% done.\n',toc(start), i/(n-1)*100)
        end
    end
    fclose(binFile);
end
% closing
for j=range
    fclose(cscFile{j});
end
fprintf(['It took ' datestr(seconds(toc(start)),'HH:MM:SS') ,' to read and write files.\n\n']);

%% Bandpass filter 600 - 6000
disp('Filtering started ...')
start2= tic;
for j=range
    [time,data,header,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc([selpath filesep 'CSC' num2str(j-1) '.ncs']);
    data_filtered = filterlfp(time, data, 600, 6000);
    write_bin_csc([selpath filesep 'CSC' num2str(j-1) '.ncs'], time,data_filtered,header,ChannelNumber,SampleFreq,NumValidSamples);
    fprintf('Filtering: %.0f seconds, %.0f%% done.\n',toc(start2), find(range==j)/length(range)*100)
end
fprintf(['It took ' datestr(seconds(toc(start2)),'HH:MM:SS') ,' to filter ', num2str(length(range)), ' csc files.\n\n']);
%% Prep for Vyash's Analysis
mkdir([selpath filesep 'B']);
ncsFiles = dir([selpath filesep '*.ncs']);
i = 0;
for j=range
    i = i+1;
    copyfile([selpath filesep 'CSC' num2str(j-1) '.ncs'], [selpath filesep 'B' filesep 'CSC_B' num2str(i) '.ncs']);
end
copyfile('VideoReport', [selpath filesep 'B' filesep 'VideoReport']);
copyfile('VideoSampling', [selpath filesep 'B' filesep 'VideoSampling']);
%% Vyash's Analysis
beep;
Vyashpath = uigetdir('C:\Users\neuropixels\Neuropixels\NPtoWinclust_Pipeline\makeParms','Select Vyash''s Code Directory');
if isa(Vyashpath,'double')
    return;
end
cd(Vyashpath)
start2= tic;
system(['ParmsGenerationPipeline.bat "' selpath filesep 'B']);
fprintf(['\nIt took ' datestr(seconds(toc(start2)),'HH:MM:SS') ,' to create the parms file.\n']);
beep;
%% Zipping
% disp('Compression started ...')
% system(['cd ' selpath '; zip All_CSCs.zip CSC*.ncs; cd ' pwd]);
%fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n\n']);
%beep;