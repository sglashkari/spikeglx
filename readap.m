%% https://djoshea.github.io/neuropixel-utils/
% Create an ImecDataset pointing at a specific
channelMapFile = '/home/shahin/Pel/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat';
imec = Neuropixel.ImecDataset('/home/shahin/Pel/07202021_1203_944_pedestal_g0_t0.exported.imec0.ap.bin', 'channelMap', channelMapFile);

meta = imec.readAPMeta();

% Reading specific time window
timeWindow = [10 30]; % in seconds
[data_partial, sampleIdx] = imec.readAP_timeWindow(timeWindow);

data_time = sampleIdx / imec.fsAP;
plot(data_time,data_partial(45,:))
                        
% Inspect the raw IMEC traces
%imec.inspectAP_timeWindow([10 20]); % 10 20 seconds into the recording