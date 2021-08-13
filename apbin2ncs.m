%% This program read a Spike GLX AP bin files to Neuralynx CSC files
% *.ap.bin ==> *.ncs
% This program is written by Shahin G Lashkari
% 2021-08-11

clc; clear
close all
range = [43:74 133:164 168:183]+1;

%%
[binaryFile,path] = uigetfile('*.ap.bin', 'Select One or More Binary Files','MultiSelect','on');
if isa(binaryFile,'double')
    return;
elseif isa(binaryFile,'char')
    binaryFile = {binaryFile};
end
selpath = uigetdir(pwd,'Select Destination Folder');
if isa(selpath,'double')
    return;
end
tic
for k = 1:length(binaryFile)
    meta = ReadMeta(binaryFile{k}, path);
    nChan = str2double(meta.nSavedChans); % 385
    nSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    n = floor(nSamp/512);
    dt = 1/str2double(meta.imSampRate); % in seconds
    tStart{k} = posixtime(datetime(strrep(meta.fileCreateTime,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone', 'America/New_York'));

    %time = 0:dt:meta.fileTimeSecs;
    
    % opening
    binFile = fopen(fullfile(path, binaryFile{k}), 'rb');
    
    if k == 1
        cscFile = cell(nChan,1); % 385 x 1
        for j=range
            FileName = ['CSC' num2str(j-1) '.ncs'];
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
        toc
    end
    
    % Body
    for i=1:n-1
        data = fread(binFile, [nChan, 512], 'int16=>double')'; % 512 x 385
        Samples = -round(data*32767 / 512); % microvolts to bits INVERTED
        TimeStamp = (tStart{k} - tStart{1} + dt + (i-1)*512*dt)*1e6; % in microseconds
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
            fwrite(cscFile{j}, Samples(:,j), 'int16');
        end
    end
    fclose(binFile);
end
% closing
for j=range
    fclose(cscFile{j});
end

toc
%%
j=46;
[time,data] = read_bin_csc(['CSC' num2str(j-1) '.ncs']);

figure(1);
plot(time,data*1e6);
xlabel('Time (sec)')
ylabel('Voltage (µV)')

%%
disp('filtering started')
for j=range
    [time,data,header,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(['CSC' num2str(j-1) '.ncs']);
    data_filtered = filterlfp(time, data, 600, 6000);
    write_bin_csc([selpath filesep 'CSC' num2str(j-1) '.ncs'], time,data_filtered,header,ChannelNumber,SampleFreq,NumValidSamples);
    fprintf('Filtering: %.0f seconds, %.0f%% done!\n',toc, find(range==j)/length(range)*100)
end
toc
%%
j=46;
[time,data] = read_bin_csc([selpath filesep 'CSC' num2str(j-1) '.ncs']);

figure(2);
plot(time,data*1e6);
xlabel('Time (sec)')
ylabel('Voltage (µV)')

%% Zipping
system(['zip ' selpath filesep 'All_CSCs_filtered_concatenated_Inverted.zip ' selpath filesep 'CSC*.ncs']);
toc
system(['rsync -tvrPh ' selpath filesep 'All_CSCs_filtered_concatenated_Inverted.zip "/run/user/1001/gvfs/smb-share:server=10.163.99.67,share=public/forShahin/"'])
toc