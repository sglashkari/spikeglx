%% https://djoshea.github.io/neuropixel-utils/
% Create an ImecDataset pointing at a specific
clear;
clc; close all;
addpath('/home/shahin/spikeglx/neuropixel-utils/+Neuropixel');
channelMapFile = '/home/shahin/spikeglx/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat';
imec = Neuropixel.ImecDataset('/home/shahin/spikeglx/07202021_1203_944_pedestal_g0_t0.exported.imec0.ap.bin', 'channelMap', channelMapFile);

meta = imec.readAPMeta();

% Reading specific time window
timeWindow = [10 10.1]; % in seconds
[data_partial, sampleIdx] = imec.readAP_timeWindow(timeWindow);

data_time = sampleIdx / imec.fsAP;  % in seconds
figure(1);
plot(data_time,data_partial(46,:));
grid on

%%
data = double(data_partial(46,:)) * 1e-6; % microVolts to Volts
size(data)
[Time,Data,Header,Samples,Timestamps] = writecsc('test.ncs', data_time, data, 0);

figure(2);
plot(Time, Data)
hold on
plot(Timestamps*1e-6, Samples(1,:),'sr')

% [Time2,Data2,Header2,Samples2, Timestamps2] = readcsc('test.ncs');
% 
% Timestamps2
% 
% figure(3);
% plot(Timestamps2*1e-6, Samples2(1,:),'sr')