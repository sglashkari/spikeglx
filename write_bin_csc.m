function write_bin_csc(filename, time,data,header,Samples,TimeStamp,ChannelNumber,SampleFreq,NumValidSamples)
%% This file is written by Shahin G. Lashkari to write Neuralynx CSC file (*.ncs)
% as binary file

if nargin < 4
    % Read a dummy header for the CSC
    headerSize = 16*1024; % 16 kbytes
    dummy = fopen('header.txt','r');
    chr = fread(dummy, [1 headerSize], '*char');
    header = splitlines(chr);
    fclose(dummy);
end


%% Write Header of CSC (size M x 1)
headerSize = 16*1024; % 16 kb
fid = fopen(filename,'w');
chr = char(join(header,newline));
chr = [chr blanks(headerSize)];
fwrite(fid, chr(1:headerSize-1), 'char');

%% Write Body of CSC
N = length(TimeStamp);
for i=1:N
    % 1 x N
    fwrite(fid, TimeStamp(1,i), 'uint64');
    fwrite(fid, ChannelNumber(1,i), 'uint32');
    fwrite(fid, SampleFreq(1,i), 'uint32');
    fwrite(fid, NumValidSamples(1,i), 'uint32');
    
    % 512 x N
    fwrite(fid, Samples(1:512,i), 'int16');
end
fclose(fid);

% %% Read header values
% InputFormat = 'yyyy/MM/dd HH:mm:ss';
% StartTimeString = extractAfter(header{8},'-TimeCreated ');
% StartTime = datetime(StartTimeString,'InputFormat',InputFormat);
% FinishTimeString = extractAfter(header{9},'-TimeClosed ');
% FinishTime = datetime(FinishTimeString,'InputFormat',InputFormat);
% 
% SamplingFrequencyString = extractAfter(header{15},'-SamplingFrequency ');
% SamplingFrequency = str2double(SamplingFrequencyString);
% 
% ADBitVoltsString = extractAfter(header{17},'-ADBitVolts ');
% ADBitVolts = str2double(ADBitVoltsString);

% %% Extract data
% data = Samples(:)* ADBitVolts; % volts
% N = size(Samples,2);
% s = 1:512:512*N; 
% sq = 1:1:512*N;
% 
% time = interp1(s,TimeStamp,sq); % check accuracy of this method
% data(isnan(time))=[];
% time(isnan(time))=[];
% time = (time * 1e-6)'; % second
% 
% %% Plot is no output is required
% if nargout == 0
%     plot(time,data);
%     fprintf('Size of data is %d.\n',length(data));
%     xlabel('Time (sec)')
%     ylabel('Voltage (V)')
%     clear time;
% end

end