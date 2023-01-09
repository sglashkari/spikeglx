clear; clc; close all
%%
[tdmsFile,tdmsPath] = uigetfile('D:\training\data.tdms','Select a tdms File to Open');
if isa(tdmsFile,'double')
    return;
end
tic

if tdmsFile == "data.tdms"
    forceFile = 'force.mat';
else
    forceFile = strrep(tdmsFile,'tdms','mat');
end

FileInfo = dir(fullfile(tdmsPath,tdmsFile));
StopTimeStamp = datetime(FileInfo.date,'InputFormat','d-MMM-y HH:mm:ss');
%%
A = TDMS_getStruct(fullfile(tdmsPath,tdmsFile));
%%
load_cell_1 = A.daq.load_cell_1.data; 
load_cell_2 = A.daq.load_cell_2.data;
load_cell_3 = A.daq.load_cell_3.data;
cam.top = A.daq.top_cam.data;
cam.side = A.daq.side_cam.data;
gap_length = A.daq.gap_length.data;

t = (1:length(load_cell_1))*5e-4; % in seconds (2kHz)
StartTimeStamp = StopTimeStamp - seconds(t(end)-t(1));
toc
%%
f1 = (load_cell_1-median(load_cell_1))/2 * 9.81; % Newton
f2 = (load_cell_2-median(load_cell_2))/2 * 9.81; % Newton
f3 = (load_cell_3-median(load_cell_3))/2 * 9.81; % Newton

%%
figure(10); clf
plot(t,cam.top>2.5,'--b'); hold on

% falling edge detection
locations = (cam.top > 2.5);
diff_locations = [0 diff(locations)];
index = (diff_locations == -1); % falling
hold on
t_fall = t(index);
plot(t_fall,zeros(size(t_fall)),'or')
fprintf('Top view camera pulses: %d \n', length(t_fall))
fprintf('On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_fall)));

t_pulse_cam_top = t_fall;
t_pulse_daq = t_pulse_cam_top; % needs to be removed in future

% side camera
plot(t,cam.side>2.5,'--g');

% falling edge detection
locations = (cam.side > 2.5);
diff_locations = [0 diff(locations)];
index = (diff_locations == -1); % falling
hold on
t_fall = t(index);
plot(t_fall,zeros(size(t_fall)),'om')
fprintf('Side view camera pulses: %d \n', length(t_fall))
fprintf('On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_fall)));
t_pulse_cam_side = t_fall;
figure(11)
plot(diff(t_pulse_cam_side),'.')
distance = gap_length;  % distance is used to be backward compatible 
save(fullfile(tdmsPath,forceFile),'t','f1','f2','f3','t_pulse_cam_top','t_pulse_cam_side','distance', 'gap_length','StartTimeStamp','StopTimeStamp');

%%
load(fullfile(tdmsPath,forceFile),'t','f1','f2','f3')
figure(1); clf
plot(t,f1,'.')
%xlim([1364 1369])
figure(2); clf
plot(t,f2,'.')
%xlim([1364 1369])
figure(3); clf
plot(t,f3,'.')
%xlim([1364 1369])

%% notch filter
fs = 1/(t(2)-t(1));
wo = 61/(fs/2);  
bw = wo;
[b,a] = iirnotch(wo,bw);
f1_filt = filtfilt(b,a,f1);
figure(1); hold on
plot(t,f1_filt)

wo = 23/(fs/2);  
bw = wo;
[b,a] = iirnotch(wo,bw);
f2_filt = filtfilt(b,a,f2);
figure(2); hold on
plot(t,f2_filt)

wo = 33/(fs/2);  
bw = wo;
[b,a] = iirnotch(wo,bw);
f3_filt = filtfilt(b,a,f3);
figure(3); hold on
plot(t,f3_filt)
save(fullfile(tdmsPath,forceFile),'f1_filt','f2_filt','f3_filt','-append');
%%
loadforce2(tdmsPath, forceFile)