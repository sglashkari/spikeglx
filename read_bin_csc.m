function [time,data,header,Samples,TimeStamp,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(filename)
%% This file is written by Shahin G. Lashkari to read  Neuralynx CSC file (*.ncs)
% as binary file
clear time data header Samples TimeStamp ChannelNumber SampleFreq NumValidSamples
if nargin == 0
     exp_directory = pwd;
     [datafile,exp_directory] = uigetfile(fullfile(exp_directory,'CSC1.ncs'), 'Select ncs File');
     filename = fullfile(exp_directory, datafile);
end

%% Read Header of CSC
headerSize = 16*1024; % 16 kbytes
fid = fopen(filename,'r');
chr = fread(fid, [1 headerSize], '*char');
header = splitlines(chr);

%% Read Body of CSC
position = ftell(fid)
i = 1;
while ~feof(fid)
    try
        % 1 x N
        TimeStamp(1,i) = fread(fid, 1, 'uint64');
        ChannelNumber(1,i) = fread(fid, 1, 'uint32');
        SampleFreq(1,i) = fread(fid, 1, 'uint32');
        NumValidSamples(1,i) = fread(fid, 1, 'uint32');
        % 512 x N
        Samples(:,i) = fread(fid, 512, 'int16');
    catch
        if ~feof(fid)
            warning('There was an error reading the file!');
        end
        break;
    end
    i = i + 1;
end
fclose(fid);

%% Read header values
InputFormat = 'yyyy/MM/dd HH:mm:ss';
StartTimeString = extractAfter(header{8},'-TimeCreated ');
StartTime = datetime(StartTimeString,'InputFormat',InputFormat);
FinishTimeString = extractAfter(header{9},'-TimeClosed ');
FinishTime = datetime(FinishTimeString,'InputFormat',InputFormat);

SamplingFrequencyString = extractAfter(header{15},'-SamplingFrequency ');
SamplingFrequency = str2double(SamplingFrequencyString);

ADBitVoltsString = extractAfter(header{17},'-ADBitVolts ');
ADBitVolts = str2double(ADBitVoltsString);

%% Extract data
if exist('Samples','var')
    data = Samples(:)* ADBitVolts; % volts
    N = size(Samples,2);
    s = 1:512:512*N;
    sq = 1:1:512*N;
    
    time = interp1(s,TimeStamp,sq); % check accuracy of this method
    data(isnan(time))=[];
    time(isnan(time))=[];
    time = (time * 1e-6)'; % second
    
    %% Plot is no output is required
    if nargout == 0
        plot(time,data);
        xlabel('Time (sec)')
        ylabel('Voltage (V)')
        clear time;
    end
else
    time=[];
    data=[];
    Samples=[];
    TimeStamp=[];
    ChannelNumber=[];
    SampleFreq=[];
    NumValidSamples=[];
end

end