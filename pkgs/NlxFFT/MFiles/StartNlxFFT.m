% clear all
if exist('NlxFFT.log', 'file')
    diary off
    delete NlxFFT.log
end
diary('NlxFFT.log')
Path = genpath(pwd);
addpath(Path);
handles = NlxFFT;
