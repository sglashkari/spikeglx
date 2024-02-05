function loadforce(tdmsPath, forceFile)
if nargin<2
    clear
    clc
    [forceFile, tdmsPath] = uigetfile('D:\training\force.mat','Select a Force File to Open');
    if isa(forceFile,'double')
        return;
    end
end

clearvars -except tdmsPath forceFile
load(fullfile(tdmsPath,forceFile),'t','f1_filt','f2_filt','f3_filt','gap_length','StartTimeStamp','StopTimeStamp')

figure(4); clf;
a(1) = subplot(2,1,1); hold on
plot(t,f1_filt)
plot(t,f3_filt)
plot(t,f2_filt)
ylabel('Load cell forces (Newtons)')
xlabel('Time (sec)')

f1p = f1_filt>2;
f2p = f2_filt>2;
f3p = f3_filt>2;
dt = 500;
shifted_f1p = [zeros(1,dt) f1p(1:end-dt)];
shifted_f2p = [zeros(1,dt) f2p(1:end-dt)];
shifted_f3p = [zeros(1,dt) f3p(1:end-dt)];

right_j = (shifted_f1p+f3p)>1;
left_j = (shifted_f3p+f1p)>1;
right_d = (shifted_f1p+f2p)>1;
left_d = (shifted_f3p+f2p)>1;

extended_right_j = ([zeros(1,dt) right_j(1:end-dt)]+right_j)>0;
extended_left_j = ([zeros(1,dt) left_j(1:end-dt)]+left_j)>0;
extended_right_d = ([zeros(1,dt) right_d(1:end-dt)]+right_d)>0;
extended_left_d = ([zeros(1,dt) left_d(1:end-dt)]+left_d)>0;

% a(2) = subplot(3,1,2); hold on
% plot(t,extended_right_j)
% plot(t,extended_left_j)
% plot(t,extended_right_d)
% plot(t,extended_left_d)

right_jump = [0 diff(extended_right_j)]>0;
left_jump = [0 diff(extended_left_j)]>0;
right_ditch = [0 diff(extended_right_d)]>0;
left_ditch = [0 diff(extended_left_d)]>0;


%% exclude double counting
% exclude double counting (jump)
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

% exclude double counting (ditch)
idx = find(right_ditch);
for i=1:nnz(right_ditch)-1
    if t(idx(i+1))-t(idx(i))<1 % 1 sec
        right_ditch(idx(i+1)) = 0;
    end
end
idx = find(left_ditch);
for i=1:nnz(left_ditch)-1
    if t(idx(i+1))-t(idx(i))<1 % 1 sec
        left_ditch(idx(i+1)) = 0;
    end
end

% exclude double counting (mistake ditch)
% time_r = sort([t(right_jump) t(right_ditch)]);
% time_l = sort([t(left_jump) t(left_ditch)]);
% t_rd = t(right_ditch);
% t_ld = t(left_ditch);
% 
% idx = find(left_ditch);
% for i=1:nnz(time_r)-1
%     ff = find(t_ld > t(idx(i)) & t_ld < t(idx(i+1)));
%      if nnz(ff) > 1
%          for j=1:nnz(ff)-1
%             plot(t_ld(idx(j)) , 0.5, 'ok')
%             left_ditch(idx(j)) = 0;
%             disp('---')
%             i
%             plot(t,left_ditch(idx(j)),'b')
%             ff(j)
%          end
%      end
% end
% 
% idx = find(right_ditch);
% for i=1:nnz(time_l)-1
%     ff = find(t_rd > t(idx(i)) & t_rd < t(idx(i+1)));
%      if nnz(ff) > 1
%          for j=1:nnz(ff)-1
%             plot(t_rd(idx(j)) , 0.5, 'ok')
%             right_ditch(idx(j)) = 0;
%             disp('---')
%             i
%             ff(j)
%          end
%      end
% end

jump = [left_jump right_jump];
ditch = [left_ditch right_ditch];
gap_length_jump = [gap_length(right_jump) gap_length(left_jump)];
gap_length_ditch = [gap_length(right_ditch) gap_length(left_ditch)];

a(2) = subplot(2,1,2); hold on
plot(t,gap_length,'k','LineWidth',2)
try
    ylim([min(gap_length)-0.5 max(gap_length)+0.5])
catch
end
ylabel('Gap length (inch)', 'FontSize', 14)
xlabel('Time (sec)', 'FontSize', 14)
grid on
plot(t,max(ylim)*right_jump,'b')
plot(t,max(ylim)*left_jump,'g')
plot(t,max(ylim)*right_ditch,'m')
plot(t,max(ylim)*left_ditch,'r')
linkaxes(a,'x')
title('Hysteresis', 'FontSize', 16)
set(gca, 'FontSize', 14)

fprintf(['\nStart date and time: ' datestr(StartTimeStamp) ,'.\n']);

fprintf(['Stop date and time: ' datestr(StopTimeStamp) ,'.\n']);
fprintf(['Total time was ' datestr(seconds(length(t)*5e-4),'HH:MM:SS') ,'.\n\n']);
fprintf('Total number of jumps: %d = %d (left) + %d (right)\n', nnz(jump), nnz(left_jump), nnz(right_jump))
if nnz(ditch) >= 5
    fprintf('Total number of ditches: %d = %d (left) + %d (right)\n', nnz(ditch), nnz(left_ditch), nnz(right_ditch))
    fprintf('Total number of passages: %d = %d (left) + %d (right) \n',  nnz(jump | ditch), nnz(left_jump+left_ditch), nnz(right_jump+right_ditch))
end
fprintf('Total jumped distance: %g inches\n', sum(gap_length_jump))
if nnz(ditch) >= 5
    fprintf('Total passed distance: %g inches\n', sum(gap_length_jump)+sum(gap_length_ditch))
end
fprintf('Distance of the jumps: [%.3f; %.3f ± %.3f; %.3f] inch\n',min(gap_length_jump), mean(gap_length_jump),std(gap_length_jump),max(gap_length_jump))
if nnz(ditch) >= 5
    fprintf('Distance of the ditches: [%.3f; %.3f ± %.3f; %.3f] inch\n',min(gap_length_ditch), mean(gap_length_ditch),std(gap_length_ditch),max(gap_length_ditch))
end
fprintf('Speed: %.f passages per hour.\n',nnz(jump | ditch)/(length(t)*5e-4)*3600)

figure(10); clf; hold on
plot(t,gap_length * 2.54,'k','LineWidth',2)
ylim([round(min(gap_length * 2.54))-1 round(max(gap_length * 2.54))+1])
ylabel('Gap length (cm)', 'FontSize', 14)
xlabel('Time (sec)', 'FontSize', 14)
grid on
set(gcf, 'Position', [100 100 1100 400]);
title(['Gap Length Change During Experiment Session, ' datestr(StartTimeStamp,'yyyy-mm-dd')], 'FontSize', 16)
set(gca, 'FontSize', 14)

% Save as a PDF
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, fullfile('E:','Rats', 'psychometric', [datestr(StartTimeStamp,'yyyy-mm-dd') '_hysteresis.pdf']), '-dpdf', '-r0')

%%
figure(5); clf
subplot(3,1,1); hold on
length_range = 3:1:30;
[hist_rj, l_r] = hist(gap_length(right_jump), length_range);
hist_rd = hist(gap_length(right_ditch), length_range);
bar(l_r, [hist_rj; hist_rd]);
title('rightward jump count')
try
    xlim([min(gap_length)-0.375 max(gap_length)+0.375])
catch
end
xlabel('Gap length (inch)')
subplot(3,1,2); hold on
[hist_lj, l_l] = hist(gap_length(left_jump), length_range);
hist_ld = hist(gap_length(left_ditch), length_range);
bar(l_l, [hist_lj; hist_ld]);
try
    xlim([min(gap_length)-0.375 max(gap_length)+0.375])
catch
end
xlabel('Gap length (inch)')
title('lefttward jump count')
subplot(3,1,3); hold on
[hist_j, l_rl] = hist(gap_length(left_jump | right_jump), length_range);
hist_d = hist(gap_length(left_ditch | right_ditch), length_range);
bar(l_rl, [hist_j; hist_d]);
try
    xlim([min(gap_length)-0.375 max(gap_length)+0.375])
catch
end
xlabel('Gap length (inch)')
title('Total passage count')
%% psychometric function (leftward vs rightward)
figure(6); clf; hold on
den = hist_rj+hist_rd;
idx = den >= 4;
plot(l_r(idx),hist_rj(idx)./(den(idx)+eps),'-o')
title('psychometric function (leftward vs rightward)')
try
    xlim([min(gap_length)-0.375 max(gap_length)+0.375])
catch
end
den = hist_lj+hist_ld;
idx = den >= 4;
plot(l_l(idx),hist_lj(idx)./(den(idx)+eps),'-o')
try
    xlim([min(gap_length)-0.375 max(gap_length)+0.375])
catch
end
den = hist_j+hist_d;
idx = den >= 4;
plot(l_rl(idx),hist_j(idx)./(den(idx)+eps),'-o')
try
    xlim([min(gap_length)-0.375 max(gap_length)+0.375])
catch
end
ylabel('Jumping percentage (%)', 'FontSize', 14)
xlabel('Gap length (inch)', 'FontSize', 14)
legend({'rightward', 'leftward', 'average'});

%% psychometric function (increase vs decrease)
figure(7); clf; hold on
% time of jump and ditch
t_rj = t(right_jump);
t_rd = t(right_ditch);
t_lj = t(left_jump);
t_ld = t(left_ditch);
t_r = [t_rj, t_rd];
t_l = [t_lj, t_ld];

rj = ones(size(t_rj));
rd = zeros(size(t_rd));
lj = ones(size(t_lj));
ld = zeros(size(t_ld));
r = [rj, rd];
l = [lj, ld];

[~, idx_r] = sort(t_r);
[~, idx_l] = sort(t_l);

jump_r = r(idx_r);
jump_l = l(idx_l);

gl_rj = gap_length(right_jump);
gl_rd = gap_length(right_ditch);
gl_lj = gap_length(left_jump);
gl_ld = gap_length(left_ditch);
gl_r = [gl_rj, gl_rd];
gl_l = [gl_lj, gl_ld];
gl_r_sorted = gl_r(idx_r);
gl_l_sorted = gl_l(idx_l);

increase_r = [0, diff(gl_r_sorted)>0];
decrease_r = [0, diff(gl_r_sorted)<0];
increase_l = [0, diff(gl_l_sorted)>0];
decrease_l = [0, diff(gl_l_sorted)<0];

% Rightward trials
subplot(2,1,1); hold on
[hist_ij, l_r] = hist(gl_r_sorted(increase_r & jump_r), length_range);
hist_id = hist(gl_r_sorted(increase_r & ~jump_r), length_range);
den = hist_ij+hist_id;
idx = den >= 2;
plot(l_r(idx),100*hist_ij(idx)./(den(idx)+eps),'rs-', 'MarkerSize', 8);  hold on

[hist_dj, l_r] = hist(gl_r_sorted(decrease_r & jump_r), length_range);
hist_dd = hist(gl_r_sorted(decrease_r & ~jump_r), length_range);
den = hist_dj+hist_dd;
idx = den >= 2;
plot(l_r(idx),100*hist_dj(idx)./(den(idx)+eps),'bo-')

x = l_r(idx);
y_increase_right = 100*hist_ij(idx)./(den(idx)+eps);
y_decrease_right = 100*hist_dj(idx)./(den(idx)+eps);
x_sigm = min(x):0.1:max(x)+1;

x100 = -10:x(1)-1;
x0 = x(end)+1:50;
x = [x100 x x0];
y_increase_right = [100*ones(size(x100)) y_increase_right zeros(size(x0))];
y_decrease_right = [100*ones(size(x100)) y_decrease_right zeros(size(x0))];
param_increase_right = sigm_fit(x, y_increase_right);
param_decrease_right = sigm_fit(x, y_decrease_right);

y_fit_increase_right = param_increase_right(1) + (param_increase_right(2)-param_increase_right(1)) ./ (1 + 10.^((param_increase_right(3)-x_sigm)*param_increase_right(4)));
y_fit_decrease_right = param_decrease_right(1) + (param_decrease_right(2)-param_decrease_right(1)) ./ (1 + 10.^((param_decrease_right(3)-x_sigm)*param_decrease_right(4)));
plot(x_sigm, y_fit_increase_right, 'r--');
plot(x_sigm, y_fit_decrease_right, 'b--');
disp('====');
disp(datestr(StartTimeStamp, 'yyyy-mm-dd'));
fprintf('Area between the data points for rightward direction: (%.2f)\n', ...
    calculate_total_area(x, y_increase_right, x, y_decrease_right)/100);
fprintf('Area between the curves for rightward direction: %.2f\n', ...
    calculate_total_area(x_sigm, y_fit_increase_right, x_sigm, y_fit_decrease_right)/100);

legend({'increase right', 'decrease right', 'fit for increase right', 'fit for decrease right'});
xlabel('Gap length (inch)', 'FontSize', 14)
ylabel('Jumping percentage (%)', 'FontSize', 14)
title('rightward (increase vs decrease)', 'FontSize', 16)
set(gcf, 'Position', [100 100 1100 700]);
h = findobj(gca, 'type', 'line');
set(h, 'LineWidth', 1.5)
set(gca, 'FontSize', 14)
ylim([0 100])


% Leftward trials
subplot(2,1,2); hold on
[hist_ij, l_l] = hist(gl_l_sorted(increase_l & jump_l), length_range);
hist_id = hist(gl_l_sorted(increase_l & ~jump_l), length_range);
den = hist_ij+hist_id;
idx = den >= 2;
plot(l_l(idx),100*hist_ij(idx)./(den(idx)+eps),'rs-', 'MarkerSize', 8); 

[hist_dj, l_l] = hist(gl_l_sorted(decrease_l & jump_l), length_range);
hist_dd = hist(gl_l_sorted(decrease_l & ~jump_l), length_range);
den = hist_dj+hist_dd;
idx = den >= 2;
plot(l_l(idx),100*hist_dj(idx)./(den(idx)+eps),'bo-')

x = l_l(idx);
y_increase_left = 100*hist_ij(idx)./(den(idx)+eps);
y_decrease_left = 100*hist_dj(idx)./(den(idx)+eps);
x_sigm = min(x):0.1:max(x)+1;

x100 = -10:x(1)-1;
x0 = x(end)+1:50;
x = [x100 x x0];
y_increase_left = [100*ones(size(x100)) y_increase_left zeros(size(x0))];
y_decrease_left = [100*ones(size(x100)) y_decrease_left zeros(size(x0))];

param_increase_left = sigm_fit(x, y_increase_left);
param_decrease_left = sigm_fit(x, y_decrease_left);
y_fit_increase_left = param_increase_left(1) + (param_increase_left(2)-param_increase_left(1)) ./ (1 + 10.^((param_increase_left(3)-x_sigm)*param_increase_left(4)));
y_fit_decrease_left = param_decrease_left(1) + (param_decrease_left(2)-param_decrease_left(1)) ./ (1 + 10.^((param_decrease_left(3)-x_sigm)*param_decrease_left(4)));
plot(x_sigm, y_fit_increase_left, 'r--');
plot(x_sigm, y_fit_decrease_left, 'b--');
fprintf('Area between the data points for leftward direction: (%.2f)\n', ...
    calculate_total_area(x, y_increase_left, x, y_decrease_left)/100);
fprintf('Area between the curves for leftward direction: %.2f\n', ...
    calculate_total_area(x_sigm, y_fit_increase_left, x_sigm, y_fit_decrease_left)/100);
disp('====');

legend({'increase left', 'decrease left', 'fit for increase left', 'fit for decrease left'});
xlabel('Gap length (inch)', 'FontSize', 14)
title('leftward (increase vs decrease)', 'FontSize', 16)
set(gcf, 'Position', [100 100 1100 700]);
h = findobj(gca, 'type', 'line');
set(h, 'LineWidth', 1.5)
ylabel('Jumping percentage (%)', 'FontSize', 14)
set(gca, 'FontSize', 14)
ylim([0 100])

%% Duplicate figure 7 to figure 8 (inch to cm)
figure(8); clf;
new_axes = copyobj(get(figure(7), 'Children'), figure(8));
set(figure(8), 'Position', [100 100 1100 700]);

% Conversion factor from inches to cm
inches_to_cm = 2.54;

% Iterate over all axes in the figure
for ax = new_axes'
    % Check if current ax is not a Legend or ColorBar
    if ~isa(ax, 'matlab.graphics.illustration.Legend') && ~isa(ax, 'matlab.graphics.illustration.ColorBar')
        % Change x and y labels to centimeters
        old_xlabels = cellfun(@str2num, get(ax, 'XTickLabel'));
        set(ax, 'XTickLabel', arrayfun(@num2str, round(old_xlabels*inches_to_cm), 'UniformOutput', false));

        % Update xlabel
        xlabel(ax, 'Gap length (cm)', 'FontSize', 14)
    end
end

% Save as a PDF
set(figure(8), 'Units', 'Inches');
pos = get(figure(8), 'Position');
set(figure(8), 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
%print(figure(8), 'hysteresis_curve.pdf', '-dpdf', '-r0')

%% Duplicate figure 4 to figure 9 (inch to cm)
figure(9); clf;
new_axes = copyobj(get(figure(4), 'Children'), figure(9));
set(figure(9), 'Position', [100 100 1100 900]);

% Conversion factor from inches to cm
inches_to_cm = 2.54;

% Iterate over all axes in the figure
for ax = new_axes'
    % Check if current ax is not a Legend or ColorBar
    if ~isa(ax, 'matlab.graphics.illustration.Legend') && ~isa(ax, 'matlab.graphics.illustration.ColorBar')
        % Change x and y labels to centimeters
        old_ylabels = cellfun(@str2num, get(ax, 'YTickLabel'));
        set(ax, 'YTickLabel', arrayfun(@num2str, round(old_ylabels*inches_to_cm), 'UniformOutput', false));
        
        % Update xlabel
        ylabel(ax, 'Gap length (cm)', 'FontSize', 16)
        hold on
    end
end

%%
if nnz(ditch) < 5
    close(6);
    close(7);
end
