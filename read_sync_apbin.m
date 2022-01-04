%% apbin2ncs reads a Spike GLX AP bin files to Neuralynx CSC files
% It also remvoes the DC offset, does common average referencing (CAR) and
% bandpass filtering (600-6000 Hz).
%
% *.ap.bin ==> 384.ncs
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

range = 385; % selection range in matlab is 1..385, while for the range for AP is 0..384

selpath = uigetdir('D:\Analysis','Select a Directory for Saving the Synchronization CSC File');
if isa(selpath,'double')
    return;
end

%% Read AP binary files and write them into csc files
chunksize = 300; % 300 samples of length 512, almost 5 seconds: 5 sec x 30Khz / 512 sample
Nlx_ADBitVolts = 0.000000036621093749999997;

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
        Samples = fread(binFile, [nChan, 512*chunksize], 'int16')'; % (512xchunksize) x 385
        Samples = de2bi(Samples(:,385),16);
        Samples = Samples(:,7)/Nlx_ADBitVolts;
        
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
                fwrite(cscFile{j}, Samples((l-1)*512+1:l*512), 'int16');
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

%% Plot
[time,data] = read_bin_csc([selpath filesep 'CSC384.ncs']);
time_np = time';
data = data/0.0012';

figure(3)
plot(time_np,data,'.')
figure(4)
plot(time_np,data<0.5,'.')

% falling edge detection
locations = (data<0.5);
diff_locations = [0;diff(locations)];
index = (diff_locations == -1);
hold on
t_fall = time_np(index);
plot(t_fall,zeros(size(t_fall)),'o')

t_np = t_fall;
fprintf('Neuropixels: On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_np)));
length(t_np)

save([selpath filesep 'np_pulse'],'time_np','data','t_np');