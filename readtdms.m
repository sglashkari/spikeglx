clear; clc; close all
[tdmsFile,tdmsPath] = uigetfile('D:\NI-Data\test\*.tdms','Tracking Data: Select a tdms File to Open');
if isa(tdmsFile,'double')
    return;
end
tdmsFile = fullfile(tdmsPath, tdmsFile);

A = TDMS_getStruct(tdmsFile);
%t = A.Untitled.Time.data;

data0 = A.Untitled.Voltage_0.data;
data1 = A.Untitled.Voltage_1.data;
data2 = A.Untitled.Voltage_2.data;
data3 = A.Untitled.Voltage_3.data;
data = [data0; data2];
t = 1:length(data);
plot(t, data,'.')
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
plot(t(index),zeros(1,sum(index)),'o')
