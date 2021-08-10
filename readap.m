%% https://djoshea.github.io/neuropixel-utils/
% Create an ImecDataset pointing at a specific

clc; clear
close all
tic
addpath([pwd filesep 'neuropixel-utils']);

channelMapFile = [pwd filesep 'neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat'];
filename = '07192021_1850_944_sleep2_g0_t0.imec0.ap.bin';
imec = Neuropixel.ImecDataset(filename, 'channelMap', channelMapFile);
meta = imec.readAPMeta();

%% Reading specific time window
timeWindow = [0 meta.fileTimeSecs]; % in seconds
[data_partial, sampleIdx] = imec.readAP_timeWindow(timeWindow);

time = sampleIdx / imec.fsAP;  % in seconds
figure(1);
plot(time,data_partial(46,:));
grid on

%%
toc
CA1 = data_partial(43:74,:);
CA3_1 = data_partial(133:164,:);
CA3_2 = data_partial(168:183,:);
toc
for i=1:32
    CA1(i,:) = filterlfp(time, CA1(i,:), 600, 6000);
    CA3_1(i,:) = filterlfp(time, CA3_1(i,:), 600, 6000);
    if i<=16
        CA3_2(i,:) = filterlfp(time, CA3_2(i,:), 600, 6000);
    end
end
toc
%%
writematrix(time,'time.csv');
writematrix(CA1,'ca1.csv');
writematrix(CA3_1,'ca31.csv');
writematrix(CA3_2,'ca32.csv');
toc
%%
% https://stackoverflow.com/questions/7618858/how-to-to-read-a-matrix-from-a-given-file
figure(2)
plot(time,CA1(1,:))