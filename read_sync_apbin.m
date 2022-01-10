function [time, data, t_fall] = read_sync_apbin(binaryFile, path)
%% read_sync_apbin reads a Spike GLX AP bin file to extract the pulse times
% used for acquiring the images from the camera
%
% *.ap.bin => 
% t_fall: time of falling signal (when the image is acquired).
%
% This program is written by Shahin G Lashkari
% Last upgrade: 2022-01-07
%
% Check:
% https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/UserManual.md#output-file-format-and-tools
%
%% pafile
if nargin < 2
    clc; clear; close all;
    [binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select a Binary File');
end
if isa(binaryFile,'double')
    return;
end
%% Read AP binary files and write them into csc files
disp('Reading synchronization data ...')
tic;
meta = ReadMeta(binaryFile, path);
nChan = str2double(meta.nSavedChans); % 385
dt = 1/str2double(meta.imSampRate); % in seconds

binFile = fopen(fullfile(path, binaryFile), 'rb');
skip = (nChan-1)*2; % int16 is 2 bytes
fseek(binFile,384*2,'bof'); % int16 is 2 bytes
data = fread(binFile, 'uint16', skip);
data = de2bi(data,16);
data = data(:,7); %7th bit
time = (0:dt:(length(data)-1)*dt)';
fclose(binFile);
fprintf('It took %.1f seconds to read synchronization data.\n\n',toc);

%%

diff_locations = [0;diff(data==1)];
t_fall = time(diff_locations == -1);

if nargout == 0
    figure(1)
    plot(time,data,'.')
    hold on
    plot(t_fall,zeros(size(t_fall)),'o')
    
    fprintf('Neuropixels: On average frames were taken every %.3f milliseconds. \n', 1000*mean(diff(t_fall)));
    fprintf('Neuropixels: Number of pulses were %d.\n',length(t_fall))
    clear time;
end
end