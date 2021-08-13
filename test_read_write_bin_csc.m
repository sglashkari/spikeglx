clc; clear;
close all;
tic
filename1 = 'CSC45.ncs';
[time,data,header,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(filename1);
% time = time(1:60*30000);
% data = (sin(2000*time)+sin(2000*time)+sin(2000*time));
plot(time,data*1e6);
xlabel('Time (sec)')
ylabel('Voltage (µV)')
ChannelNumber(1:10)
SampleFreq(1:10);
%%
filename2 = 'test.ncs';
%data_filtered = filterlfp(time, data, 600, 6000);
write_bin_csc(filename2, time,data,header,ChannelNumber,SampleFreq,NumValidSamples);
[time,data,header,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(filename2);
ChannelNumber(1:10)
SampleFreq(1:10);
toc
figure(1)
hold on
plot(time,data*1e6,'o');
xlabel('Time (sec)')
ylabel('Voltage (µV)')
%xlim([1205 1205.2])