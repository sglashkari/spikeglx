clc; clear
close all
addpath([pwd filesep 'pkgs/releaseDec2015/binaries']); % Neuralynx packages for Linux/Mac
addpath([pwd filesep 'neuropixel-utils']);
addpath([pwd filesep 'neuropixel-utils/map_files']);
%%
%[Time1,Data1,Header1,Samples1,Timestamps1,Data_bits1] = readcsc('CSC1.ncs',[1000 1100]*1e6);
% figure(1)
% plot(Time1,Data1)
% figure(2)
% plot(Timestamps1*1e-6, Samples1(1,:),'sr')

%%

channelMapFile = [pwd filesep 'neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat'];

filename = '07192021_1850_944_sleep2_g0_t0.imec0.ap.bin';

imec = Neuropixel.ImecDataset(filename, 'channelMap', channelMapFile);

meta = imec.readAPMeta();

% Reading specific time window
timeWindow = [0 meta.fileTimeSecs]; % in seconds
[data_partial, sampleIdx] = imec.readAP_timeWindow(timeWindow);

data_time = sampleIdx / imec.fsAP;  % in seconds
figure(1);
plot(data_time,data_partial(46,:));
grid on

%% Write to CSC
for i = 1:meta.nSavedChans
    time = data_time'; % in seconds
    data = double(data_partial(i,:))' * 1e-6; % microVolts to Volts
    
    % time = Time1;
    % data = Data1;
    
    FileName = ['CSC' num2str(i-1) '.ncs'];
    AppendToFileFlag = 0; % Delete the file if exists
    ExportMode = 1; % Export All
    ExportModeVector = 1; % (Extract All): The vector value is ignored.
    %FieldSelectionFlags = [1 1 1 1 1 1]; % Timestamps, Channel Numbers, Sample Frequency, Number of Valid Samples, Samples, Header
    FieldSelectionFlags = [1 0 0 0 1 0];
    %
    
    %   Timestamps: A 1xN integer vector of timestamps. This must be in ascending order.
    N = floor(length(time)/512);
    time = time(1:512*N)';
    Timestamps = time(1:512:end) * 1e6; % microseconds
    
    % Number Of Records In Matlab Arrays
    NumRecs = N;
    
    %   Samples: A 512xN integer matrix of the data points. These values are in AD counts.
    data_bits = round(data(1:512*N) / 0.000000036621093749999997); % volts to bits
    Samples = reshape(data_bits,512,N);
    
    FieldSelection = [1 0 0 0 1 0];
    
    system(['rm -f ' FileName]);
    
    Mat2NlxCSC(FileName, AppendToFileFlag, ExportMode, ExportModeVector, NumRecs, FieldSelectionFlags, Timestamps, Samples);
    
end

figure(2)
plot(time,data_bits)
hold on
plot(Timestamps*1e-6, Samples(1,:),'sr')

%% Read From CSC
FieldSelectionFlags = [1 0 0 0 1]; % Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples
HeaderExtractionFlag = 0;
ExtractionMode = 1;
ExtractionModeVector = [];

[Timestamps, Samples] = Nlx2MatCSC_v3(FileName, FieldSelectionFlags, HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);

%%
Data = Samples(:)* 0.000000036621093749999997; % volts
N = size(Samples,2);
s = 1:512:512*N; 
sq = 1:1:512*N;

Time = interp1(s,Timestamps,sq); % check accuracy of this method

Data(isnan(Time))=[];
Time(isnan(Time))=[];

%Time = Timestamps(1) + (1:length(Data)) * 1e6 / 30000; % micrseconds
Time = (Time * 1e-6)'; % seconds

% [Time3,Data3,Header3,Samples3,Timestamps3,Data_bits3] = readcsc('test.ncs');
figure(3)
plot(Time, Data*1e6)
ylabel('Voltage (µV)')
xlabel('Time (sec)')

% figure(4)
% plot(Timestamps*1e-6, Samples(1,:),'sr')