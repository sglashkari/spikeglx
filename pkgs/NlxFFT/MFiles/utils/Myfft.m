function [S, F] = Myfft(s, Fs, varargin)
% [S, F] = Myfft(s, Fs, |h, REAL, dB, NFFT|)
% |h, REAL, dB, N| are optional
% fft S(F) = s(Fs)
% h is either a figure handle or an axes handle to plot to.
% if REAL is true, default, then return the abs of the first 1/2 of fft
% dB(true/false) is only relevant if REAL = true
% NFFT = 2^N, default is 2^nextpow2(length of s)

s = squeeze(s);
if ndims(s) > 2
    error('first arg must be either matrix or squeezable to matrix')
end
if isvector(s)
    s = s(:);
end
    
if length(varargin) < 4 || isempty(varargin{4})
    NFFT = 2^nextpow2(size(s,1));
else
    NFFT = varargin{4};
    if mod(NFFT,2) 
        %Force NFFT even
        NFFT = NFFT+1;
    end
end


if length(varargin) < 3 || isempty(varargin{3})
    dB = true;
else
    dB = varargin{3};
end

if length(varargin) < 2 || isempty(varargin{2})
    REAL = true;
else
    REAL = varargin{2};
end

if length(varargin) < 1 || isempty(varargin{1})
    h = false;
else
    h = varargin{1};
end


%% see http://www.ni.com/white-paper/4278/en
% for list of windows advantages, disadvantages and corresponding scaling
% factors. Note also the amplitude of an fft will always be a bit off due
% to leakage.
HammingScaleFactor = 0.54;
F = linspace(0, Fs, NFFT);
Window = repmat(hamm(size(s,1))/HammingScaleFactor, 1, size(s, 2));

s = Window.*s;
S = 2*fft(s, NFFT)/length(s);

if REAL
    S = abs(S(1:NFFT/2,:));
    F = F(1:NFFT/2);
    if dB
        S = 20*log10(S);
    end
end
if h && ishandle(h)
    if ~isreal(S)
        warning('In general, if plotting, set 4th arg, REAL, to true or leave blank.')
    end
    switch get(h, 'type')
        case 'axes'
            plot(h, F, S);
        case 'figure'
            figure(h);
            plot(F, S);
        otherwise
            warning('First arg must be [ ], axes or figure handle');
    end
    
    xlabel('Frequency')
    ylabel('Volt')
end