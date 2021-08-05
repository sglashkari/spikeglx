%% https://djoshea.github.io/neuropixel-utils/
% Create an ImecDataset pointing at a specific
addpath('/home/shahin/spikeglx/neuropixel-utils/+Neuropixel');
channelMapFile = '/home/shahin/spikeglx/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat';
imec = Neuropixel.ImecDataset('/home/shahin/spikeglx/07202021_1203_944_pedestal_g0_t0.exported.imec0.ap.bin', 'channelMap', channelMapFile);

meta = imec.readAPMeta();

% Reading specific time window
timeWindow = [10 30]; % in seconds
[data_partial, sampleIdx] = imec.readAP_timeWindow(timeWindow);

data_time = sampleIdx / imec.fsAP;
plot(data_time,data_partial(45,:));
