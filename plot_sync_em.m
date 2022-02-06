% Test and plot sync_em
% SGL 2022-02-05
clear; 
clc; close all
%% Selecting the appropriate files
[binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select a Binary File');
if isa(binaryFile,'double')
    return;
end

[csvFile,csvPath] = uigetfile('full-tracking.csv','Tracking Data: Select a CSV File to Open');
if isa(csvFile,'double')
    return;
end

[forceFile,forcePath] = uigetfile('D:\NI-Data\force.mat','Select a force.mat File to Open');
if isa(forceFile,'double')
    return;
end

csvFile = fullfile(csvPath, 'full-tracking.csv');
forceFile = fullfile(forcePath,'force.mat');
%% Neuropixels
[time,data,t_pulse_np] = read_sync_apbin(binaryFile, path);

figure(1)
plot(time,data,'.')
hold on
plot(t_pulse_np,zeros(size(t_pulse_np)),'o')

fprintf('Neuropixels: On average frames were taken every %.6g milliseconds. \n', 1000*mean(diff(t_pulse_np)));
fprintf('Neuropixels: Number of pulses were %d.\n',length(t_pulse_np));
fprintf('Neuropixels: length of pulses recording %4.10g seconds.\n\n',t_pulse_np(end)-t_pulse_np(1));

%% Camera
disp('Reading tracking data ...')
start = tic;
T = readtable(csvFile);
t_cam = T.t;

fprintf('Camera: On average frames were taken every %.6g milliseconds. \n', 1000*mean(diff(t_cam)));
fprintf('Camera: Number of frames taken were %d.\n',length(t_cam))
fprintf('Camera: length of recording %4.12g seconds.\n\n',t_cam(end)-t_cam(1));

%% DAQ (load cells)
load(forceFile,'t','t_pulse_daq')

fprintf('DAQ: On average frames were taken every %.6g milliseconds. \n', 1000*mean(diff(t_pulse_daq)));
fprintf('DAQ: Number of frames taken were %d.\n',length(t_pulse_daq))
fprintf('DAQ: length of recording %4.12g seconds.\n',t_pulse_daq(end)-t_pulse_daq(1));

%% Comparison
N = min(length(t_pulse_np),length(t_cam));

figure(2); clf;
plot(1:length(t_pulse_np),t_pulse_np-t_pulse_np(1),'.',1:N,t_cam-t_cam(1),'.')
legend('np','camera')
xlabel('t np')

p = polyfit(t_pulse_np(1:N),t_pulse_np(1:N)-t_cam(1:N),1);
t_diff = polyval(p,t_pulse_np(1:N));
%t_diff = 0;

figure(3); clf;
plot(t_cam(1:N),t_pulse_np(1:N)-t_cam(1:N),'.')
title('diff np and camera')

figure(4); clf;
plot(2:length(t_pulse_np),diff(t_pulse_np(1:end)),'.',2:N,diff(t_cam(1:N)),'.')
legend('np','camera')

%% Neuropixels: plot changes of diff

figure(5); clf;
% edge detection
singal_np = [0; diff(t_pulse_np)];
singal_np = movmedian(singal_np,10);
locations_np = double(singal_np>0.032);
diff_locations_np = [0; diff(locations_np)];
index = (diff_locations_np == -1 | diff_locations_np == 1); % falling or rising
hold on
t_fall_np = t_pulse_np(index);
plot(t_fall_np,zeros(size(t_fall_np)),'o')

plot(t_pulse_np,locations_np)
length(t_fall_np)

%% Camera: plot changes of diff

figure(5); hold on;
t_pulse_cam = t_cam;
% edge detection
singal_cam = [0; diff(t_pulse_cam)];
singal_cam = movmedian(singal_cam,10);
locations_cam = double(singal_cam>0.032);
diff_locations_cam = [0; diff(locations_cam)];
index = (diff_locations_cam == -1 | diff_locations_cam == 1); % falling or rising
hold on
t_fall_cam = t_pulse_cam(index);
plot(t_fall_cam,zeros(size(t_fall_cam)),'o')

plot(t_pulse_cam,locations_cam)
length(t_fall_cam)

%% interpolation
t_cam_in_np = sync_em(t_pulse_np,t_pulse_cam);

figure(6); clf;
plot(t_pulse_np,zeros(size(t_pulse_np)),'.b'); hold on
plot(t_cam_in_np,zeros(size(t_cam_in_np)),'or');
plot(t_fall_np,zeros(size(t_fall_np)),'pk');

%%
fprintf('Camera (corrected): length of recording %4.12g seconds.\n\n',t_cam_in_np(end)-t_cam_in_np(1));
