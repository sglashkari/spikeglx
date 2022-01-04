clear; clc; close all
%%
[tdmsFile,tdmsPath] = uigetfile('D:\NI-Data\test\data.tdms','Tracking Data: Select a tdms File to Open');
if isa(tdmsFile,'double')
    return;
end
tic
tdmsFile = fullfile(tdmsPath, tdmsFile);
A = TDMS_getStruct(tdmsFile);
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
data = [data0; data1-1; data2];
sync = [data3; data4];
figure;
plot(t, sync)

figure;
plot(t, data,'-.')
% hold on
% data_filtered = filterlfp(t,data,0.1,30);
% plot(t, data_filtered)
figure
plot(t,data3>2.5,'.')


% falling edge detection
locations = (data3 <= 2.5);
diff_locations = [0 diff(locations)];
index = (diff_locations == 1);
hold on
t_fall = t(index);
plot(t_fall,zeros(size(t_fall)),'o')
fprintf('On average frames were taken every %.6f milliseconds. \n', 1000*mean(diff(t_fall)));
length(t_fall)

plot(t, data3,'-.'); hold on; plot(t, data,'-.');
%xlim([1262.2 1362.5])
%%
try
    data5 = A.Untitled.Untitled_5.data;
%    data6 = A.Untitled.Untitled_6.data;
%    data7 = A.Untitled.Untitled_7.data;
%    rat = A.Untitled.Untitled_8.data;
catch
end

ax1 = subplot(4,1,1);
plot(t,data)
ax2 = subplot(4,1,2);
plot(t,data3)
ax3 = subplot(4,1,3);
plot(t,data4)
ax4 = subplot(4,1,4);
plot(t,data5)
linkaxes([ax1 ax2 ax3 ax4],'x')
