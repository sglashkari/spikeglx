function [Time,Data,Header,Samples] = readcsc(ncs_filename, TimeRange)
if nargin == 0
     exp_directory = pwd;
     [datafile,exp_directory] = uigetfile(fullfile(exp_directory,'*.ncs'), 'Select ncs File');
     ncs_filename = fullfile(exp_directory, datafile);
end
FieldSelectionFlags = [1 1 1 1 1]; % Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples
HeaderExtractionFlag = 1;

if nargin < 2
    ExtractionMode = 1;
    ExtractionModeVector = [];
else
    ExtractionMode = 4;
    ExtractionModeVector = TimeRange;
end

if ispc
    addpath('pkgs/MatlabImportExport_v6.0.0'); % Neuralynx packages for Windows
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
        Samples, Header] = Nlx2MatCSC(ncs_filename, FieldSelectionFlags,...
        HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
else
    addpath('pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
        Samples, Header] = Nlx2MatCSC_v3(ncs_filename, FieldSelectionFlags,...
        HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
end

SamplingFrequency = 30000;
ADBitVolts = 0.000000036621093749999997;

Data = Samples(:)*ADBitVolts; % volts
N = size(Samples,2);
s = 1:512:512*N; 
sq = 1:1:512*N;

Time = interp1(s,Timestamps,sq); % check accuracy of this method

Data(isnan(Time))=[];
Time(isnan(Time))=[];

Time = Timestamps(1) + (1:length(Data)) * 1e6 / SamplingFrequency; % micrseconds
Time = (Time * 1e-6)'; % seconds

if nargout == 0
    plot(Time,Data);
    fprintf('Size of data is %d.\n',length(Data));
    clear Time;
end
    
end