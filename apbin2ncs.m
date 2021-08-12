%% This program read a Spike GLX AP bin files to Neuralynx CSC files
% *.ap.bin ==> *.ncs
% This program is written by Shahin G Lashkari
% 2021-08-11

clc; clear
close all
tic
range = [43:74 133:164 168:183]+1;

%%
[binaryFile,path] = uigetfile('*.ap.bin', 'Select One or More Binary Files','MultiSelect','on');
if isa(binaryFile,'double')
    return;
elseif isa(binaryFile,'char')
    binaryFile = {binaryFile};
end

for k = 1:length(binaryFile)
    meta = ReadMeta(binaryFile{k}, path);
    nChan = str2double(meta.nSavedChans); % 385
    nSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    n = floor(nSamp/512);
    dt = 1/str2double(meta.imSampRate); % in seconds
    %time = 0:dt:meta.fileTimeSecs;
    
    % opening
    binFile = fopen(fullfile(path, binaryFile{k}), 'rb');
    
    if k == 1
        cscFile = cell(nChan,1); % 385 x 1
        for j=range
            FileName = ['CSC' num2str(j-1) '.ncs'];
            cscFile{j} = fopen(FileName,'w');
        end
        
        % Header
        headerSize = 16*1024; % 16 kb
        dummy = fopen('header.txt','r');
        hdr = fread(dummy, [1 headerSize], '*char');
        fclose(dummy);
        for j=range
            fwrite(cscFile{j}, hdr, 'char'); % Write Header of CSCs
        end
        toc
    end
    
    % Body
    for i=1:n-1
        Samples = fread(binFile, [nChan, 512], 'int16')'; % 512 x 385
        TimeStamp = (dt + (i-1)*512*dt)*1e6; % in microseconds
        SampleFreq=str2double(meta.imSampRate);
        NumValidSamples=512;
        for j=range
            ChannelNumber=j-1;
            % Write Body of CSC
            fwrite(cscFile{j}, TimeStamp, 'uint64');
            fwrite(cscFile{j}, ChannelNumber, 'uint32');
            fwrite(cscFile{j}, SampleFreq, 'uint32');
            fwrite(cscFile{j}, NumValidSamples, 'uint32');
            % 512 x N
            fwrite(cscFile{j}, Samples(1:512,j), 'int16');
        end
    end
    fclose(binFile);
end
% closing
for j=range
    fclose(cscFile{j});
end


toc

j=46;
[time,data,header,Samples,TimeStamp,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(['CSC' num2str(j-1) '.ncs']);

figure;
plot(time,data);
xlabel('Time (sec)')
ylabel('Voltage (V)')