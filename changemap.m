clear
load('neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat')
range = 43:74;
chanMap = chanMap(range);
chanMap0ind = chanMap0ind(range);
connected = connected(range);
shankInd = shankInd(range);
xcoords = xcoords(range);
ycoords = ycoords(range);
save('neuropixel-utils/map_files/CA1.mat')