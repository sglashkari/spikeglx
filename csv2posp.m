%% pos_p_gen generates Pos.p from a given tracking data
% full tracking.csv ==> pos.p & pos.p.ascii
% It creates pos.p files for all the subfolders
% 2022-01-07
% Author Shahin G Lashkari
clear;
clc; close all
pixels_width = 640;
pixels_height = 480;

%% Selecting the appropriate files
[binaryFile,path] = uigetfile('D:\NeuralData\*.ap.bin', 'Select a Binary File');
if isa(binaryFile,'double')
    return;
end

[csvFile,csvPath] = uigetfile('full-tracking.csv','Tracking Data: Select a CSV File to Open');
if isa(csvFile,'double')
    return;
end

selparentpath = uigetdir('D:\Analysis','Select the Main Experiment Directory for Saving Pos.p Files');
if isa(selparentpath,'double')
    return;
end

Vyashpath = 'C:\Users\neuropixels\Neuropixels\NPtoWinclust_Pipeline\makeParms';
Vyashpath = uigetdir(Vyashpath,'Select Vyash''s Code Directory');
if isa(Vyashpath,'double')
    return;
end

answer=questdlg('Is the tracking time in the CSV already synchronized with Neuropixels time?', 'Warning!', 'Yes', 'No', 'No');

thisFilePath = pwd;
paramspath = {dir(fullfile(selparentpath,'**','*.parms')).folder}';
[paramspath, idx] = unique(paramspath);
paramsfilename = {dir(fullfile(selparentpath,'**','*.parms')).name}';
paramsfilename = paramsfilename(idx);

%% Read tracking data from file
% The table should have t, x, y, (hd) as headers
disp('Reading tracking data ...')
start = tic;
csvFile = fullfile(csvPath, csvFile);
T = readtable(csvFile);
try
    idx = T.x>=0; % condition for successful tracking (this need to be adjusted depending on the tracker)
    Tf = T(idx,:); % filtered table
    x = Tf.x;
    y = Tf.y;
catch
    disp('The CSV file should have at least 3 columns with t, x, y & hd (hd is optional) as column headers.');
    return;
end
try
    angle = Tf.hd;
    hd_exists = true;
catch
    disp('Head directions are not detected!');
    angle = -99*ones(size(x));
    hd_exists = false;
end

%% Time
if isequal(answer,'Yes')
    disp('Tracking data are already in Neuropixels time!')
    ts = Tf.t * 1e6; % in microseconds
else
    [~,~,t_pulse_np] = read_sync_apbin(binaryFile, path);
    if length(t_pulse_np)==height(T)
        fprintf('Number of pulses in Neuropixels and camera were the same: %d.\n',length(t_pulse_np))
    else
        warning('Number of pulses do not match between Neuropixels and camera.');
        fprintf('Neuropixels: Number of pulses were %d.\n',length(t_pulse_np))
        fprintf('Camera: Number of frames taken were %d.\n',height(T))
    end
    disp('Tracking data is syncronized with Neuropixels time!')
    ts = (Tf.t - Tf.t(1) + t_pulse_np(1)) * 1e6; % in microseconds
end

%% Calculations
box_x = [min(x) max(x)];
box_y = [min(y) max(y)];

scale_x = pixels_width / diff(box_x);
scale_y = pixels_height / diff(box_y);
scale = min([scale_x scale_y])*0.95;

box_center = mean([box_x; box_y],2);

pos_x = (x- box_center(1)) * scale + pixels_width/2;
pos_y = (y- box_center(2)) * scale + pixels_height/2;

figure; plot(pos_x, pixels_height - pos_y,'.')
axis equal; axis([0 pixels_width 0 pixels_height])

figure; plot(ts*1e-6,pos_x,'.b');
xlabel('Time (sec)'); ylabel('X position');ylim([0 pixels_width])

if hd_exists
    figure; plot(ts*1e-6,angle,'.b');
    xlabel('Time (sec)'); ylabel('Head direction (deg)')
end

%% Writing Pos.p
disp('Starting ...')
meta = ReadMeta(binaryFile, path);
start_time = datetime(strrep(meta.fileCreateTime,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone', 'America/New_York');

header = {"%%BEGINHEADER";...
    "%Tracking data from DeepLabCut:";...
    strcat("%Date of the Experiment: ",datestr(start_time,"mmmm dd, yyyy"));...
    strcat("%Time of the Experiment: ",datestr(start_time,"HH:MM AM"));...
    "%File Type: Binary";...
    "%Format:";...
    "%--------------------------";...
    "%|  Timestamp   (Double)  |";...
    "%|------------------------|";...
    "%|  X coord     (Single)  |";...
    "%|------------------------|";...
    "%|  Y coord     (Single)  |";...
    "%|------------------------|";...
    "%|  Direction   (Single)  |";...
    "%--------------------------";...
    "%";...
    "%Use following record type to read in:";...
    "%Public Type PosRecord";...
    "%    TS As Double";...
    "%    X As Single";...
    "%    Y As Single";...
    "%    Dir As Single";...
    "%End Type";...
    "%%ENDHEADER"};

for l = 1:length(paramspath)
    cd(thisFilePath);
    fprintf('\n%d from %d: \n\n',l,length(paramspath));
    %% Writing header of the files
    
    posppath = fullfile(paramspath{l},'..');
    pos_p_file_name = fullfile(posppath,'Pos.p');
    pos_p_ascii_file_name = fullfile(posppath,'Pos.p.ascii');
    
%     if isfile(pos_p_file_name)
%         beep;
%         answer=questdlg('The Pos.p file already exists! Do you want to overwrite it?', 'Warning!', 'Yes', 'No', 'Yes');
%         if ~isequal(answer,'Yes')
%             disp('Pos.p was not overwritten!')
%             return
%         end
%     end
    
    file = fopen(pos_p_file_name,'w');
    file_ascii = fopen(pos_p_ascii_file_name,'w');
    format_spec = '%s\r\n';
    
    for row = 1:size(header,1)
        fprintf(file,format_spec,header{row,:});
        fprintf(file_ascii,format_spec,header{row,:});
    end
    
    %% Writing data to the file in binary and ascii format
    positions = [ts, pos_x, pos_y, angle];
    format_spec_ascii = '%.0f, %.4f, %.4f, %.4f\r\n';
    
    for row = 1:size(positions,1)
        fwrite(file,ts(row),'float64');
        fwrite(file,positions(row,2:4),'float32');
        fprintf(file_ascii,format_spec_ascii,positions(row,:));
    end
    
    fclose(file);
    disp([pos_p_file_name ' created!'])
    fclose(file_ascii);
    disp([pos_p_ascii_file_name ' created!'])
    
    %% Vyash's Analysis
    cd(Vyashpath);
    pause(1);
    start2= tic;
    % writing parms_B_withPos.Ntt.parms and parms_B_withPos.waveforms
    system(['python addDLCPosToParms.py ' fullfile(paramspath{l},paramsfilename{l}) ' ' pos_p_ascii_file_name]);
    parmsname = extractBefore(paramsfilename{l},'.Ntt.parms');
    movefile(fullfile(paramspath{l},[parmsname '.waveforms']),fullfile(paramspath{l},[parmsname '_withPos.waveforms']));
    delete(fullfile(paramspath{l},paramsfilename{l}));
end
%%
cd(thisFilePath);
fprintf(['It totally took ' datestr(seconds(toc(start)),'HH:MM:SS') ,'.\n']);