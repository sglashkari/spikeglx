function loadforce(tdmsPath, forceFile)
if nargin<2
    [forceFile, tdmsPath] = uigetfile('D:\training\force.mat','Select a Force File to Open');
    if isa(forceFile,'double')
        return;
    end
end

clearvars -except tdmsPath forceFile
load(fullfile(tdmsPath,forceFile),'t','f1_filt','f2_filt','f3_filt','gap_length','StartTimeStamp','StopTimeStamp')

figure(4); clf; clear a
a(1) = subplot(2,1,1);
plot(t,[f1_filt;f3_filt])
ylabel('Load cell forces (Newtons)')
xlabel('Time (sec)')

f1p = f1_filt>2;
f2p = f2_filt>2;
f3p = f3_filt>2;
dt = 500;
shifted_f1p = [zeros(1,dt) f1p(1:end-dt)];
shifted_f2p = [zeros(1,dt) f2p(1:end-dt)];
shifted_f3p = [zeros(1,dt) f3p(1:end-dt)];

right = (shifted_f1p+(f3_filt>2))>1;
left = (shifted_f3p+(f1_filt>2))>1;
extended_right = ([zeros(1,dt) right(1:end-dt)]+right)>0;
extended_left = ([zeros(1,dt) left(1:end-dt)]+left)>0;
right_jump = [0 diff(extended_right)]>0;
left_jump = [0 diff(extended_left)]>0;
% exclude double counting
idx = find(right_jump);
for i=1:nnz(right_jump)-1
    if t(idx(i+1))-t(idx(i))<1 % 1 sec
        right_jump(idx(i+1)) = 0;
    end
end
idx = find(left_jump);
for i=1:nnz(left_jump)-1
    if t(idx(i+1))-t(idx(i))<1 % 1 sec
        left_jump(idx(i+1)) = 0;
    end
end

a(2) = subplot(2,1,2); hold on
plot(t,gap_length,'k')
try
    ylim([min(gap_length)-0.5 max(gap_length)+0.5])
catch
end
ylabel('Gap length (inch)')
xlabel('Time (sec)')
grid on
plot(t,max(ylim)*right_jump,'r')
plot(t,max(ylim)*left_jump,'g')
linkaxes(a,'x')
fprintf(['\nStart date and time: ' datestr(StartTimeStamp) ,'.\n']);
fprintf(['Stop date and time: ' datestr(StopTimeStamp) ,'.\n']);
fprintf(['Total time was ' datestr(seconds(t(end)),'HH:MM:SS') ,'.\n\n']);
fprintf('Total number of jumps: %d = %d (right) + %d (left)\n', nnz(left_jump+right_jump), nnz(left_jump), nnz(right_jump))
fprintf('Total jumped distance: %g inches\n', sum(gap_length(right_jump))+sum(gap_length(left_jump)))
fprintf('Distance of the right jumps: [%.3f; %.3f ± %.3f; %.3f] inch\n',min(gap_length(right_jump)), mean(gap_length(right_jump)),std(gap_length(right_jump)),max(gap_length(right_jump)))
fprintf('Distance of the left jumps: [%.3f; %.3f ± %.3f; %.3f] inch\n',min(gap_length(left_jump)), mean(gap_length(left_jump)),std(gap_length(left_jump)),max(gap_length(left_jump)))
fprintf('Speed: %.f jumps per hour.\n',nnz(left_jump+right_jump)/t(end)*3600)

%%
figure(5); clf
subplot(2,1,1)
length_range = 3:0.25:30;
histogram(gap_length(right_jump), length_range)
title('rightward jump count')
try
    xlim([min(min(gap_length(right_jump)),min(gap_length(left_jump)))-0.125 max(max(gap_length(right_jump)),max(gap_length(left_jump)))+0.375])
catch
end
xlabel('Gap length (inch)')
subplot(2,1,2)
histogram(gap_length(left_jump), length_range)
try
    xlim([min(min(gap_length(right_jump)),min(gap_length(left_jump)))-0.125 max(max(gap_length(right_jump)),max(gap_length(left_jump)))+0.375])
catch
end
xlabel('Gap length (inch)')
title('leftward jump count')
