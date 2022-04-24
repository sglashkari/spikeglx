%%WINCLUST2MAT extracts data from the output of winclust (cluster mazes) and
% tracking.dat (tracking data) and save it in data.mat file in the
% experiment directory.
%
%   See also LAPDETECTOR, ANALYZEDATA.
%
% Date 2022-01-16 (originally 2021-01-31)
% Author Shahin G Lashkari
%
clc; clear; close all;
%% Selecting the appropriate files
[binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select a Binary File');
if isa(binaryFile,'double')
    return;
end

[csvFile,csvPath] = uigetfile('Z:\Rat980\full-tracking.csv','CSV File: Select a Tracking Data to Open');
if isa(csvFile,'double')
    return;
end

% [sideTimeFile,sideTimePath] = uigetfile(fullfile(csvPath,'side_times.csv'),'Side Camera Times: Select a CSV File to Open');
% if isa(sideTimeFile,'double')
%     return;
% end

csvFile = fullfile(csvPath, 'full-tracking.csv');
sideTimeFile = fullfile(csvPath, 'side_times.csv');
%dlcFilename = fullfile(csvPath, 'dlc.csv');

[forceFile,forcePath] = uigetfile('D:\NI-Data\force.mat','Select a mat File to Open');
if isa(forceFile,'double')
    return;
end

selparentpath = uigetdir('D:\Analysis','Select the Main Experiment Directory');
if isa(selparentpath,'double')
    return;
end

%% Experiment information

% finding start and finish time
meta = ReadMeta(binaryFile, path);
exp.start = datetime(strrep(meta.fileCreateTime,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone', 'America/New_York');
exp.finish = exp.start + str2double(meta.fileTimeSecs)/24/3600;
exp.date = datetime(extractBefore(meta.fileCreateTime,'T'),'InputFormat','yyyy-MM-dd','TimeZone', 'America/New_York');

% finding rat number and day number
if exp.date < datetime(2021,12,31,'TimeZone','America/New_York') &&  exp.date > datetime(2021,11,10,'TimeZone','America/New_York')
    exp.name = datestr(exp.date);
    exp.rat_no = 980;
    exp.day = day(exp.date);
else
    error('specify the rat number!');
end

%% Neural Data
disp('Reading clusters data ...')
start = tic;
listing = dir(fullfile(selparentpath,'**','cl-maze*.*'));
folders = {listing.folder}';
names = {listing.name}';
absolue_paths = fullfile(folders, names);
N = length(listing); % number of clusters

sh_no = str2double(extractBetween(folders,'Shank',filesep));
region = char(convertCharsToStrings(extractBetween(folders,[selparentpath filesep],'-Shank'))); % like CA1, CA3, PPC
maze_no = str2double(extractBetween(names,'maze','.'));
cluster_no = str2double(extractAfter(names,'.'));
cluster_no(isnan(cluster_no))=0;    % cluster 0

A = cellfun(@(x) importdata(x,',',13), absolue_paths, 'UniformOutput', false);
fprintf(['It took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);

%% Camera data
% Read tracking data from file
% The table should have t, x, y, (hd) as headers
T = readtable(csvFile);
idx1 = T.x>=0; % condition for successful tracking => idx1: index for single-marker tracker
%(this need to be adjusted depending on the tracker)

%figure(1); clf; plot(T.frame_no(idx),T.x(idx),'.r');

% Exclude bad tracking data points
bad_frames = [];
switch exp.date
    case '10-Dec-2021'
        bad_frames = [45202:45289 45967:46066 47688:47755]; % [inital final)
end
idx1(bad_frames)=0;
%hold on; plot(T.frame_no(idx),T.x(idx),'.k'); pause;

% position
% 2020-03 5 ft = 1524 mm = 480 pixels (each pixel = 3.175 mm)
% 2020-10 3 ft = 914 mm = 840 pixels = norm([296 372]-[1136 348],2) >> each pixel ~ 1.1 mm
% 2021-12 3 ft = 914 mm = 662 pixels = norm([1366 229.3]-[704 206.156]) >> each pixel ~ 1.4 mm

if exp.rat_no >= 980
    ppcm = norm([1366 229.3]-[704 206.156])/91.4; % pixels per cm
else
    ppcm = norm([296 372]-[1136 348],2)/91.4; % pixels per cm
end

pos.x = T.x(idx1) / ppcm; % cm
pos.y = T.y(idx1) / ppcm; % cm

% time
[~,~,t_pulse_np] = read_sync_apbin(binaryFile, path);
if length(t_pulse_np)==height(T)
    fprintf('Number of pulses in Neuropixels and camera were the same: %d.\n',length(t_pulse_np));
    pos.time = t_pulse_np; % in seconds
else
    warning('Number of pulses do not match between Neuropixels and camera.');
    fprintf('Neuropixels: Number of pulses were %d.\n',length(t_pulse_np));
    fprintf('Camera: Number of frames taken were %d.\n',height(T));
    pos.time = sync_em(t_pulse_np,T.t); % in seconds
end
pos.t=pos.time(idx1);
pos.frame = T.frame_no(idx1); % frame number starts from 0
pos.len = T.ditch_length(idx1);

% velocity
pos.vx = gradient(pos.x)./gradient(pos.t); % Vx in cm/sec
pos.vy = gradient(pos.y)./gradient(pos.t); % Vy in cm/sec
pos.s = vecnorm([pos.vx pos.vy]')'; % speed in cm/sec
%pos.hd = atan2d(pos.vy,pos.vx); % estimation of hd based
%pos.ax = gradient(pos.vx)./gradient(pos.t); % ax in cm/sec

ax1 = subplot(2,1,1);
plot(pos.t,pos.x,'.')
ylabel('x')
ax2 = subplot(2,1,2);
plot(pos.t,pos.vx,'.')
ylabel('vx')
linkaxes([ax1 ax2],'x')

%% Balazs's tracking (position and quaternion)
close all;
pos.p = [T.pos_1 T.pos_2 T.pos_3];
pos.q = [T.rot_4 T.rot_1 T.rot_2 T.rot_3];
pos.success = T.success;
pos.p = pos.p(idx1,:);
pos.q = pos.q(idx1,:);
pos.success = pos.success(idx1);
idx2 = pos.success==1; % idx2: index for the crown tracker

[yaw, pitch, roll] = quat2angle(pos.q); % idx1: index for single-marker tracker
pos.yaw = rad2deg(yaw);
pos.pitch = rad2deg(pitch);
pos.roll = rad2deg(roll);

% correction for the miscalibration of the table
coefficients = polyfit(pos.p(idx2,1),pos.p(idx2,3), 1);
phi = atan(coefficients(1))*180/pi;
R = [cosd(phi) 0 sind(phi); 0 1 0; -sind(phi) 0 cosd(phi)];
figure
subplot(2,1,1); plot(pos.p(idx2,1),pos.p(idx2,3),'.', 'MarkerSize',0.2);
pos.p = (100*R*pos.p')';  % 100 for m to cm
subplot(2,1,2); plot(pos.p(idx2,1),pos.p(idx2,3),'.', 'MarkerSize',0.2);

% exclude the tracking data when the crown is occluded as NaN
figure
plot3(pos.p(:,1),pos.p(:,2),pos.p(:,3),'.'); hold on
if strcmp('21-Dec-2021',exp.name)
    idx3 = pos.p(:,2)< -22.5 | pos.p(:,2) > 25 | pos.p(:,3)> 24 | pos.p(:,3)< -35;
    idx5 = [0; abs(diff(pos.p(:,1))) > 10 | abs(diff(pos.p(:,2))) > 5];
else
    idx3 = zeros(size(idx2));
    idx5 = idx3;
end
idx4 = ~idx2|idx3|idx5; % need to be excluded

figure(2)
pos.p = interp1(pos.t(~idx4),pos.p(~idx4,:),pos.t);
plot3(pos.p(idx4,1),pos.p(idx4,2),pos.p(idx4,3),'or'); hold on

pos.p(idx4,1:3)=nan;
figure(2); plot3(pos.p(:,1),pos.p(:,2),pos.p(:,3),'.g'); hold on

%% DAQ data (load cell, side view camera pulses)
try % for temporary fix before updating force.mat files
    load(fullfile(forcePath,'force.mat'),'t','f1','f2','f3','t_pulse_cam_top','t_pulse_cam_side','f1_filt','f2_filt','f3_filt');
    t_pulse_cam_top = reshape(t_pulse_cam_top,[],1);
    t_pulse_cam_side = reshape(t_pulse_cam_side,[],1);
catch
    load(fullfile(forcePath,'force.mat'),'t','f1','f2','f3','t_pulse_daq','f1_filt','f2_filt','f3_filt');
    t_pulse_cam_top = t_pulse_daq;
end
% time
if length(t_pulse_np)==length(t_pulse_cam_top)
    fprintf('Number of pulses in Neuropixels and NI DAQ were the same: %d.\n',length(t_pulse_np));
    daq.t = interp1(t_pulse_cam_top,t_pulse_np,t,'linear','extrap'); % in seconds
else
    warning('Number of pulses do not match between Neuropixels and NI DAQ.');
    fprintf('Neuropixels: Number of pulses were %d.\n',length(t_pulse_np));
    fprintf('NI DAQ: Number of pulses were %d.\n',length(t_pulse_cam_top));
    t_pulse_daq_in_np = sync_em(t_pulse_np,t_pulse_cam_top); % in seconds
    daq.t = interp1(t_pulse_cam_top,t_pulse_daq_in_np,t,'linear','extrap'); % in seconds
end
daq.loadcell = [f1;f2;f3];
daq.filt.loadcell = [f1_filt;f2_filt;f3_filt];

%% Camera (side view)
T_side = readmatrix(sideTimeFile);
t_cam_side = round(T_side(:,5),3); % seconds (with milliseconds precision)
try
    t_side_cam_in_daq = sync_em(t_pulse_cam_side,t_cam_side);
    pos.side.t = interp1(t, daq.t, t_side_cam_in_daq,'linear','extrap'); % side view camera time in neuropixel time
    pos.side.frame = 0:length(pos.side.t)-1;
catch
end

%% Deep lab cut (side view)
% try
%     dlc = readtable(dlcFilename);
%     dlc.z = dlc.y;
%     dlc.t = pos.side.t(ismember(dlc.frame,pos.side.frame));
%     figure;
%     plot(dlc.t,dlc.z);
%     title('Elevation')
%     xlim('Time (sec)')
% catch
%     dlc = [];
% end
%% spike data
cluster(N).name ='';
for idx2 = 1:N
    cluster(idx2).name = [region(idx2,:) '_shank' num2str(sh_no(idx2)) '_cluster' num2str(cluster_no(idx2))];
    cluster(idx2).region = region(idx2,:);
    cluster(idx2).sh = sh_no(idx2);
    cluster(idx2).m = maze_no(idx2);
    cluster(idx2).cl = cluster_no(idx2);
    cluster(idx2).no = idx2;
    % ineterpolation for time (excluding times that the rat is occluded)
    cluster(idx2).t = (A{idx2}.data(:,18))*1e-6; % sec
    % interpolation for position
    cluster(idx2).x = interp1(pos.t, pos.x, cluster(idx2).t);
    cluster(idx2).y = interp1(pos.t, pos.y, cluster(idx2).t);
    cluster(idx2).p = interp1(pos.t, pos.p, cluster(idx2).t);
    % interpolation for velocity
    cluster(idx2).vx = interp1(pos.t, pos.vx, cluster(idx2).t);
    cluster(idx2).vy = interp1(pos.t, pos.vy, cluster(idx2).t);
    cluster(idx2).s = vecnorm([cluster(idx2).vx cluster(idx2).vy]')';
end

%% Saving
mat_filename = fullfile(selparentpath,'data.mat');
save(mat_filename,'pos','cluster','exp','ppcm', 'daq');
disp(['File ' mat_filename ' has been created!'])
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);