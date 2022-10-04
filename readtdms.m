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

A = TDMS_getStruct(fullfile(tdmsPath,tdmsFile));
%%
data0 = A.daq.load_cell_1.data;
data1 = A.daq.load_cell_2.data;
data2 = A.daq.load_cell_3.data;
cam.top = A.daq.top_cam.data;
cam.side = A.daq.side_cam.data;
distance = A.daq.gap_length.data;
t = (1:length(data0))*5e-4; % in seconds (2kHz)
StartTimeStamp = StopTimeStamp - seconds(t(end)-t(1));
toc
%%
f1 = data0/2 * 9.81; % Newton
f2 = (data1-median(data1))/2 * 9.81; % Newton
f3 = data2/2 * 9.81; % Newton

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
try
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

    save(fullfile(tdmsPath,forceFile),'t','f1','f2','f3','t_pulse_daq','t_pulse_cam_top','t_pulse_cam_side','distance'); % t_pulse_daq needs to be removed in future
catch
    save(fullfile(tdmsPath,forceFile),'t','f1','f2','f3','t_pulse_daq'); % needs to be removed in future
end
%xlim([1262.2 1362.5])
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
load(fullfile(tdmsPath,forceFile),'t','f1_filt','f2_filt','f3_filt','distance')
figure(4); clf; clear a
a(1) = subplot(2,1,1);
plot(t,[f1_filt;f3_filt])
ylabel('Load cell forces (Newtons)')
xlabel('Time (sec)')

f1p = f1_filt>2;
f2p = f2_filt>2;
f3p = f3_filt>2;
dt = 500;
shifted_f1p = [zeros(1,dt) f1p(1:end-dt)];
shifted_f2p = [zeros(1,dt) f2p(1:end-dt)];
shifted_f3p = [zeros(1,dt) f3p(1:end-dt)];

right = (shifted_f1p+(f3_filt>2))>1;
left = (shifted_f3p+(f1_filt>2))>1;
extended_right = ([zeros(1,dt) right(1:end-dt)]+right)>0;
extended_left = ([zeros(1,dt) left(1:end-dt)]+left)>0;
right_jump = [0 diff(extended_right)]>0;
left_jump = [0 diff(extended_left)]>0;
% exclude double counting
idx = find(right_jump);
for i=1:nnz(right_jump)-1
    if t(idx(i+1))-t(idx(i))<1 % 1 sec
        right_jump(idx(i+1)) = 0;
    end
end
idx = find(left_jump);
for i=1:nnz(left_jump)-1
    if t(idx(i+1))-t(idx(i))<1 % 1 sec
        left_jump(idx(i+1)) = 0;
    end
end

a(2) = subplot(2,1,2); hold on
plot(t,distance,'k')
try
    ylim([min(min(distance(right_jump)),min(distance(left_jump)))-0.5 max(max(distance(right_jump)),max(distance(left_jump)))+0.5])
catch
end
ylabel('Gap length (inch)')
xlabel('Time (sec)')
grid on
plot(t,max(ylim)*right_jump,'r')
plot(t,max(ylim)*left_jump,'g')
linkaxes(a,'x')
fprintf(['\nStart date and time: ' datestr(StartTimeStamp) ,'.\n']);
fprintf(['Stop date and time: ' datestr(StopTimeStamp) ,'.\n']);
fprintf(['Total time was ' datestr(seconds(t(end)),'HH:MM:SS') ,'.\n\n']);
fprintf('Total number of jumps: %d = %d (right) + %d (left)\n', nnz(left_jump+right_jump), nnz(left_jump), nnz(right_jump))
fprintf('Total jumped distance: %g inches\n', sum(distance(right_jump))+sum(distance(left_jump)))
fprintf('Distance of the right jumps: [%.3f; %.3f ± %.3f; %.3f] inch\n',min(distance(right_jump)), mean(distance(right_jump)),std(distance(right_jump)),max(distance(right_jump)))
fprintf('Distance of the left jumps: [%.3f; %.3f ± %.3f; %.3f] inch\n',min(distance(left_jump)), mean(distance(left_jump)),std(distance(left_jump)),max(distance(left_jump)))
fprintf('Speed: %.f jumps per hour.\n',nnz(left_jump+right_jump)/t(end)*3600)

%%
figure(5); clf
subplot(2,1,1)
length_range = 3:0.25:30;
histogram(distance(right_jump), length_range)
title('rightward jump count')
try
    xlim([min(min(distance(right_jump)),min(distance(left_jump)))-0.125 max(max(distance(right_jump)),max(distance(left_jump)))+0.375])
catch
end
xlabel('Gap length (inch)')
subplot(2,1,2)
histogram(distance(left_jump), length_range)
try
    xlim([min(min(distance(right_jump)),min(distance(left_jump)))-0.125 max(max(distance(right_jump)),max(distance(left_jump)))+0.375])
catch
end
xlabel('Gap length (inch)')
title('leftward jump count')
