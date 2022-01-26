%%WINCLUST2MAT extracts data from the output of winclust (cluster mazes) and
% tracking.dat (tracking data) and save it in data.mat file in the
% experiment directory.
%
%   See also LAPDETECTOR, ANALYZEDATA.
%
%   SGL 2022-01-16 (originally 2021-01-31)
%

clc; clear; close all;
%% Selecting the appropriate files
[binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select a Binary File');
if isa(binaryFile,'double')
    return;
end

[csvFile,csvPath] = uigetfile('full-tracking.csv','CSV File: Select a Tracking Data to Open');
if isa(csvFile,'double')
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
region = char(convertCharsToStrings(extractBetween(folders,[selparentpath filesep],'-Shank'))); % like CA1, CA3, ...
maze_no = str2double(extractBetween(names,'maze','.'));
cluster_no = str2double(extractAfter(names,'.'));
cluster_no(isnan(cluster_no))=0;    % cluster 0

A = cellfun(@(x) importdata(x,',',13), absolue_paths, 'UniformOutput', false);
fprintf(['It took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);

%% Position data
% Read tracking data from file
% The table should have t, x, y, (hd) as headers
csvFile = fullfile(csvPath, csvFile);
T = readtable(csvFile);
idx = T.x>=0; % condition for successful tracking (this need to be adjusted depending on the tracker)
Tf = T(idx,:); % filtered table

% position
% 2020-03 5 ft = 1524 mm = 480 pixels (each pixel = 3.175 mm)
% 2020-10 3 ft = 914 mm = 840 pixels = norm([296 372]-[1136 348],2) >> each pixel ~ 1.1 mm
% 2021-12 3 ft = 914 mm = 662 pixels = norm([1366 229.3]-[704 206.156]) >> each pixel ~ 1.4 mm
ppcm = norm([1366 229.3]-[704 206.156])/91.4; % pixels per cm
pos.x = Tf.x / ppcm; % cm
pos.y = Tf.y / ppcm; % cm

% time
[~,~,t_pulse_np] = read_sync_apbin(binaryFile, path);
if length(t_pulse_np)==height(T)
    fprintf('Number of pulses in Neuropixels and camera were the same: %d.\n',length(t_pulse_np));
    offset = t_pulse_np - T.t;
    pos.t = t_pulse_np; % in seconds
else
    warning('Number of pulses do not match between Neuropixels and camera.');
    fprintf('Neuropixels: Number of pulses were %d.\n',length(t_pulse_np));
    fprintf('Camera: Number of frames taken were %d.\n',height(T));
    offset = t_pulse_np(1)- T.t(1);
    pos.t = T.t + offset; % in seconds
end
pos.t=pos.t(idx);
pos.frame = Tf.frame_no; % frame number starts from 0

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

%% spike data
cluster(N).name ='';
for index = 1:N
    cluster(index).name = [region(index,:) '_shank' num2str(sh_no(index)) '_cluster' num2str(cluster_no(index))];
    cluster(index).region = region(index,:);
    cluster(index).sh = sh_no(index);
    cluster(index).tt = cluster(index).sh;  % sh for shank same as tt for tetrode
    cluster(index).m = maze_no(index);
    cluster(index).cl = cluster_no(index);
    cluster(index).no = index;
    %cluster(index).ti = str2double(A{index}.textdata{12})*1e-6; % sec
    %cluster(index).tf = str2double(A{index}.textdata{13})*1e-6; % esc
    cluster(index).t = (A{index}.data(:,18))*1e-6; % sec
    % interpolation for position
    cluster(index).x = interp1(pos.t, pos.x, cluster(index).t);
    cluster(index).y = interp1(pos.t, pos.y, cluster(index).t);
    % interpolation for velocity
    cluster(index).vx = interp1(pos.t, pos.vx, cluster(index).t);
    cluster(index).vy = interp1(pos.t, pos.vy, cluster(index).t);
    cluster(index).s = vecnorm([cluster(index).vx cluster(index).vy]')';
end

%% Saving
mat_filename = fullfile(selparentpath,'data.mat');
save(mat_filename,'pos','cluster','exp','ppcm', 'offset');
disp(['File ' mat_filename ' has been created!'])
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);