function t_cam_in_np = sync_em(t_pulse_np,t_pulse_cam)
%
% This program synchronizes two clocks:
%   e.g. the Neuropixels clock and the camera clock.
%
%   See also PLOT_SYNC_EM, PLOT_SYNC_SIDE.
%
%   Date 2022-02-05
%
%   Author Shahin G Lashkari
%

%% Camera: change detection for variable frame rate
singal_cam = [0; diff(t_pulse_cam)];
singal_cam = movmedian(singal_cam,10);
locations_cam = double(singal_cam>0.032);
diff_locations_cam = [0; diff(locations_cam)];
index = (diff_locations_cam == -1 | diff_locations_cam == 1); % falling or rising
t_fall_cam = t_pulse_cam(index);

%% Neuropixels: change detection for variable frame rate
singal_np = [0; diff(t_pulse_np)];
singal_np = movmedian(singal_np,10);
locations_np = double(singal_np>0.032);
diff_locations_np = [0; diff(locations_np)];
index = (diff_locations_np == -1 | diff_locations_np == 1); % falling or rising
t_fall_np = t_pulse_np(index);

%% interpolation
if length(t_fall_np) == length(t_fall_cam)
    t_cam_in_np = interp1(t_fall_cam, t_fall_np, t_pulse_cam,'linear','extrap');
    fprintf('Sync was successful!\n\n');
else
    t_cam_in_np = interp1(t_fall_cam, t_fall_np(1:length(t_fall_cam)), t_pulse_cam,'linear','extrap');
    warning('Different number of pulses; Clocks may not be synced!')
end

end