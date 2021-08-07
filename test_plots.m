%% Read From CSC
i = 320;
FileName = ['CSC' num2str(i) '.ncs'];

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