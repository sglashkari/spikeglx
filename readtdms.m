clear; clc; close all
%%
[tdmsFile,tdmsPath] = uigetfile('D:\NI-Data\test\data.tdms','Select a tdms File to Open');
if isa(tdmsFile,'double')
    return;
end
%%
tic
A = TDMS_getStruct(fullfile(tdmsPath,'data.tdms'));
toc

try % before Dec 13, 2021
    data0 = A.Untitled.Voltage_0.data;
    data1 = A.Untitled.Voltage_1.data;
    data2 = A.Untitled.Voltage_2.data;
    data3 = A.Untitled.Voltage_3.data;
    t = (1:length(data0))*1e-3; % in seconds (1kHz)
catch
    try
        data0 = A.Untitled.PXI1Slot4_ai0.data;
        data1 = A.Untitled.PXI1Slot4_ai1.data;
        data2 = A.Untitled.PXI1Slot4_ai2.data;
        data3 = A.Untitled.PXI1Slot4_ai3.data;
        t = (1:length(data0))*1e-4; % in seconds
    catch
        data0 = A.Untitled.Dev1_ai0.data;
        data1 = A.Untitled.Dev1_ai1.data;
        data2 = A.Untitled.Dev1_ai2.data;
        data3 = A.Untitled.Dev1_ai3.data;
        t = (1:length(data0))*5e-4; % in seconds (2kHz)
        try
            data4 = A.Untitled.Dev1_ai4.data;
        catch
            fprintf('No AI 4!\n');
            t = (1:length(data0))*1e-4; % in seconds (10kHz)
        end
    end
end

%%
figure(10)
plot(t,data3>2.5,'.')

% falling edge detection
locations = (data3 > 2.5);
diff_locations = [0 diff(locations)];
index = (diff_locations == -1); % falling
hold on
t_fall = t(index);
plot(t_fall,zeros(size(t_fall)),'o')
fprintf('On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_fall)));
length(t_fall)

%% 
f1 = data0/2 * 9.81; % Newton
f2 = (data1-median(data1))/2 * 9.81; % Newton
f3 = data2/2 * 9.81; % Newton
t_pulse_daq = t_fall;
save(fullfile(tdmsPath,'force.mat'),'t','f1','f2','f3','t_pulse_daq');
toc
%xlim([1262.2 1362.5])
%%
tic
load(fullfile(tdmsPath,'force.mat'),'t','f1','f2','f3','t_pulse_daq')
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