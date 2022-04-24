clear; clc; close all
%%
[tdmsFile,tdmsPath] = uigetfile('D:\NI-Data\test\data.tdms','Select a tdms File to Open');
if isa(tdmsFile,'double')
    return;
end
tic
A = TDMS_getStruct(fullfile(tdmsPath,'data.tdms'));
toc
%%
try % before Dec 13, 2021
    data0 = A.Untitled.Voltage_0.data;
    data1 = A.Untitled.Voltage_1.data;
    data2 = A.Untitled.Voltage_2.data;
    cam.top = A.Untitled.Voltage_3.data;
    t = (1:length(data0))*1e-3; % in seconds (1kHz)
catch
    try
        data0 = A.Untitled.PXI1Slot4_ai0.data;
        data1 = A.Untitled.PXI1Slot4_ai1.data;
        data2 = A.Untitled.PXI1Slot4_ai2.data;
        cam.top = A.Untitled.PXI1Slot4_ai3.data;
        t = (1:length(data0))*1e-4; % in seconds
    catch
        data0 = A.Untitled.Dev1_ai0.data;
        data1 = A.Untitled.Dev1_ai1.data;
        data2 = A.Untitled.Dev1_ai2.data;
        cam.top = A.Untitled.Dev1_ai3.data;
        t = (1:length(data0))*5e-4; % in seconds (2kHz)
        try
            cam.side = A.Untitled.Dev1_ai4.data;
        catch
            fprintf('No AI 4!\n');
            t = (1:length(data0))*1e-4; % in seconds (10kHz)
        end
    end
end

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
fprintf('On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_fall)));
length(t_fall)
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
    fprintf('On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_fall)));
    length(t_fall)
    t_pulse_cam_side = t_fall;
    figure(11)
    plot(diff(t_pulse_cam_side),'.')

    save(fullfile(tdmsPath,'force.mat'),'t','f1','f2','f3','t_pulse_daq','t_pulse_cam_top','t_pulse_cam_side'); % t_pulse_daq needs to be removed in future
catch
    save(fullfile(tdmsPath,'force.mat'),'t','f1','f2','f3','t_pulse_daq'); % needs to be removed in future
end
toc
%xlim([1262.2 1362.5])
%%
tic
load(fullfile(tdmsPath,'force.mat'),'t','f1','f2','f3')
toc
figure(1); clf
plot(t,f1,'.')
%xlim([1364 1369])
figure(2); clf
plot(t,f2,'.')
%xlim([1364 1369])
figure(3); clf
plot(t,f3,'.')
%xlim([1364 1369])
figure(4); clf
plot(t,[f1;f2;f3])
toc

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

%%
save(fullfile(tdmsPath,'force.mat'),'f1_filt','f2_filt','f3_filt','-append');