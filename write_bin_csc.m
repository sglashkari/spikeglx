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

%% Read header values
ADBitVoltsString = extractAfter(header{17},'-ADBitVolts ');
ADBitVolts = str2double(ADBitVoltsString);

%% Extract data

%   Timestamps: A 1xN integer vector of timestamps. This must be in ascending order.
N = floor(length(time)/512);
time = time(1:512*N)';
TimeStamp = time(1:512:end) * 1e6; % microseconds

%   Samples: A 512xN integer matrix of the data points. These values are in AD counts.
data_bits = round(data(1:512*N) / ADBitVolts); % volts to bits
Samples = reshape(data_bits,512,N);

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

end