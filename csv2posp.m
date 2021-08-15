%pos_p_gen generates Pos.p from a given tracking data
% *.csv ==> pos.p & pos.p.ascii
% August 13, 2021
% Author Shahin G Lashkari
clear; 
clc; close all

pixels_width = 640;
pixels_height = 480;
frame_rate = 30; % 30 fps
offset = 0;
%% Read data from file
[csvFile,path] = uigetfile('*.csv','Select CSV File to Open');
if isa(csvFile,'double')
    return;
end

csvFile = fullfile(path, csvFile);
A = readmatrix(csvFile,'Range','A:O');

p = A(:,4);
position_flag = p>0.95;

t = A(position_flag,1)/frame_rate; % in seconds
x = A(position_flag,2);
y = A(position_flag,3);

hd = A(position_flag,14);
angle = filterhd(hd);

%% Calculations
ts = t * 1e6 + offset; % in microseconds

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

figure; plot(t,pos_x,'.b'); 
xlabel('Time (sec)'); ylabel('X position');ylim([0 pixels_width])
figure; plot(t,angle,'.b'); 
xlabel('Time (sec)'); ylabel('Head direction (deg)')

%% Writing header of the files
start_time = '';
header = {"%%BEGINHEADER";...
       "%Stadium-style rig for jumping rats: Created from 2D head tracking data";...
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

pos_p_file_name = fullfile(path,'Pos.p');
pos_p_ascii_file_name = fullfile(path,'Pos.p.ascii');

if isfile(pos_p_file_name)
    answer=questdlg('The Pos.p file already exists! Do you want to overwrite it?', 'Warning!', 'Yes', 'No', 'Yes');
    if ~isequal(answer,'Yes')
        disp('Pos.p was not overwritten!')
        return
    end
end

file = fopen(pos_p_file_name,'w');
file_ascii = fopen(pos_p_ascii_file_name,'w');
format_spec = '%s\r\n';

for row = 1:size(header,1)
    fprintf(file,format_spec,header{row,:});
    fprintf(file_ascii,format_spec,header{row,:});
end

%% Writing data to the file in binary and ascii format
positions = [ts, pos_x, pos_y, angle];
format_spec_ascii = '%.2f, %.4f, %.4f, %.4f\r\n';

for row = 1:size(positions,1)
    fwrite(file,ts(row),'float64');
    fwrite(file,positions(row,2:4),'float32');
    fprintf(file_ascii,format_spec_ascii,positions(row,:));
end

fclose(file);
disp('Pos.p created!')
fclose(file_ascii);
disp('Pos.p.ascii created!')
