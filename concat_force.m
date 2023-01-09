clear
clc
[forceFile, tdmsPath] = uigetfile('D:\training\force.mat','Select Force Files to Concatenate', 'MultiSelect', 'on');
if isa(forceFile,'double')
    return;
end

clearvars -except tdmsPath forceFile
load(fullfile(tdmsPath,forceFile{1}),'t','f1_filt','f2_filt','f3_filt','gap_length','StartTimeStamp','StopTimeStamp')
combined.t = t;
combined.f1_filt = f1_filt;
combined.f2_filt = f2_filt;
combined.f3_filt = f3_filt;
combined.gap_length = gap_length;
combined.StartTimeStamp = StartTimeStamp;
combined.StopTimeStamp = StopTimeStamp;

for i=2:length(forceFile)
    load(fullfile(tdmsPath,forceFile{i}),'t','f1_filt','f2_filt','f3_filt','gap_length','StartTimeStamp','StopTimeStamp')
    combined.t = [combined.t t+seconds(StartTimeStamp-combined.StopTimeStamp)];
    combined.f1_filt = [combined.f1_filt f1_filt];
    combined.f2_filt = [combined.f2_filt f2_filt];
    combined.f3_filt = [combined.f3_filt f3_filt];
    combined.gap_length = [combined.gap_length gap_length];
    combined.StopTimeStamp = StopTimeStamp;
end


t = combined.t;
f1_filt = combined.f1_filt;
f2_filt = combined.f2_filt;
f3_filt = combined.f3_filt;
gap_length = combined.gap_length;
StartTimeStamp = combined.StartTimeStamp;
StopTimeStamp = combined.StopTimeStamp;

save(fullfile(tdmsPath,['concat_', forceFile{1}]),'t','f1_filt','f2_filt','f3_filt', 'gap_length','StartTimeStamp','StopTimeStamp');
loadforce2(tdmsPath, ['concat_', forceFile{1}])