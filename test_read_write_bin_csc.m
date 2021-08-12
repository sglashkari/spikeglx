clc; clear;
close all;
filename = 'CSC1.ncs';
[time,data,header,Samples,TimeStamp,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(filename);

plot(time,data);
fprintf('Size of data is %d.\n',length(data));
xlabel('Time (sec)')
ylabel('Voltage (V)')

%%
filename = 'test.ncs';
write_bin_csc(filename, time,data,header,Samples,TimeStamp,ChannelNumber,SampleFreq,NumValidSamples);
[time,data,header,Samples,TimeStamp,ChannelNumber,SampleFreq,NumValidSamples] = read_bin_csc(filename);

figure(2)
plot(time,data);
fprintf('Size of data is %d.\n',length(data));
xlabel('Time (sec)')
ylabel('Voltage (V)')