function writecsc(ncs_filename, time, data, ChNo)
%WRITECSC 
%

% % old
% Filename = [extractBefore(Filename,'_00') '.ntt'];
% AppendToFileFlag = 0; % Delete file if exists
% ExportMode = 1; % Export All
% ExportModeVector = 1; % (Extract All): The vector value is ignored.
% NumRecs = length(Timestamps); % Number of records in arrays
% FieldSelectionFlags = [1 1 1 1 1 1]; % Timestamps, Spike Channel Numbers, Cell Numbers, Features, Samples, Header
% 
% Mat2NlxTT( Filename, AppendToFileFlag, ExportMode, ExportModeVector, NumRecs, ...
%     FieldSelectionFlags, Timestamps, ScNumbers, CellNumbers, ...
%     Features, Samples, Header);

% new
FileName = ncs_filename;
AppendToFileFlag = 0; % Delete the file if exists
ExportMode = 1; % Export All
ExportModeVector = 1; % (Extract All): The vector value is ignored.
FieldSelectionFlags = [1 1 1 1 1 1]; % Timestamps, Channel Numbers, Sample Frequency, Number of Valid Samples, Samples, Header

%   Timestamps: A 1xN integer vector of timestamps. This must be in ascending order.
N = floor(length(time)/512);
Timestamps = time(1:512:end);

%   ChannelNumbers: A 1xN integer vector of channel numbers.
ChannelNumbers = ChNo * ones(1,N); % between 0 to 31

%   SampleFrequencies: A 1xN integer vector of sample frequencies.
SampleFrequencies = 30000 * ones(1,N);

%   NumberOfValidSamples: A 1xN integer vector of the number of valid samples in the
%                         corresponding item in the Sample output variable.
NumberOfValidSamples = 512 * ones(1,N);

%   Samples: A 512xN integer matrix of the data points. These values are in AD counts.
Samples = reshape(data(1:512*N),512,N);






%% {'-ADBitVolts 0.000000036621093749999997'} volts to bits = 1/ ADBitVolts;
%% Time = (Time * 1e-6)'; % seconds








%   Header: A Mx1 string vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.
Header = {'This file has created from a spike GLX file'};

Mat2NlxCSC(FileName, AppendToFileFlag, ExportMode, ExportModeVector, ...
    FieldSelectionFlags, Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
                Samples, Header);

end

