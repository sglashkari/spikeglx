function writecsc(ncs_filename, time, data, ChNo) 
%WRITECSC 
%
%   data in Volts
%   time in seconds
%   Channel Number 0 .. 31
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
Timestamps = time(1:512:end) * 1e6; % microseconds

% Number Of Records In Matlab Arrays
NumRecs = N;

%   ChannelNumbers: A 1xN integer vector of channel numbers.
ChannelNumbers = ChNo * ones(1,N); % between 0 to 31

%   SampleFrequencies: A 1xN integer vector of sample frequencies.
SampleFrequencies = 30000 * ones(1,N);

%   NumberOfValidSamples: A 1xN integer vector of the number of valid samples in the
%                         corresponding item in the Sample output variable.
NumberOfValidSamples = 512 * ones(1,N);

%   Samples: A 512xN integer matrix of the data points. These values are in AD counts.
data = round(data / 0.000000036621093749999997); % volts to bits
Samples = reshape(data(1:512*N),512,N);

%   Header: A Mx1 string vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.
Header = {'######## This file has created from a spike GLX file';
        '-FileType NCS'                                                  ;
        '-FileVersion 3.4'                                               ;
        '-FileUUID e86b5a32-eca0-47d7-9105-28cfbc302e55'                 ;
        '-SessionUUID cab85976-6efd-436e-88f6-f32b6af5cd44'              ;
        '-ProbeName '                                                    ;
        '-OriginalFileName "C:\CheetahData\2020-11-22_18-32-37\CSC4.ncs"';
        '-TimeCreated 2020/11/22 18:33:37'                          ;     
        '-TimeClosed 2020/11/22 19:44:50'                          ;      
        ''                                                        ;
        '-RecordSize 1044'                                               ;
        '-ApplicationName Cheetah "6.4.0 Development"'                   ;
        '-AcquisitionSystem AcqSystem1 Cube2'                            ;
        '-ReferenceChannel "External Hardware Reference"'                ;
        '-SamplingFrequency 30000'                ;                       
        '-ADMaxValue 32767'                      ;                        
        '-ADBitVolts 0.000000036621093749999997'};

    disp('hello')
if ispc
    addpath('pkgs/MatlabImportExport_v6.0.0'); % Neuralynx packages for Windows
    Mat2NlxCSC(FileName, AppendToFileFlag, ExportMode, ExportModeVector, ...
                FieldSelectionFlags, Timestamps, ChannelNumbers, SampleFrequencies, ...
                NumberOfValidSamples,Samples, Header);
else
    addpath('pkgs/releaseDec2015/binaries'); % Neuralynx packages for Linux/Mac
    Mat2NlxCSC(FileName, AppendToFileFlag, ExportMode, ExportModeVector, NumRecs, ...
                FieldSelectionFlags, Timestamps, ChannelNumbers, SampleFrequencies, ...
                NumberOfValidSamples, Samples, Header);            
end

end