function [rhythm, phase, mag]= filterlfp(t, lfp, w, wch)
% FILTERLFP filters out theta from the raw LFP signal
% 
% (I)   filterlfp(time, signal, w)
%
% Two methods:
%
% (1) w = {'lfp'; 'ap'; 'delta'; 'theta'; 'beta'; 'gamma'}
% (2) w = [wcl wch] e.g [6 9] for theta
%       wcl : lower cutoff frequency
%       wch : higher cutoff frequency
%
% (II)  filterlfp(time, signal, wcl, wch)
%
%   Date 2023-01-03
%
%   Author Shahin G Lashkari
%
if nargin < 2
    [t, lfp] = read_bin_csc;
end

Ts = mean(diff(t));
SamplingFrequency = round(1/Ts);

if nargin < 3
    w = 'theta';
elseif nargin == 4
    w = [w wch];
end

if isstring(w) || ischar(w)
    switch w
        case 'lfp'
            w = [1 400];
        case 'ap'
            w = [600 6000];     
        case 'delta'
            w = [0.5 3.5];
        case 'theta'
            w = [6 9];
        case 'beta'
            w = [10 20];
        case 'gamma'
            w = [30 50];
        otherwise
            w = [6 9];
            warning('Type is not recognized, but theta is chosen!')
    end
end
wn = w/(SamplingFrequency/2);
order = 2;

[b,a] = butter(order, wn);
rhythm = filtfilt(b,a,double(lfp));

if nargout == 0
    close all;
    plot(t, lfp, 'Color', uint8([230 230 230]));
    hold on
    plot(t, rhythm,'r');
    clear rhythm;
end
if nargout > 1
    z = hilbert(rhythm);
    phase = rad2deg(angle(z)); 
end
if nargout > 2
    mag = abs(z);
end
    
end