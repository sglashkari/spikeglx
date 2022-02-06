function hd_filtered = filterhd(hd)
if nargin == 0
    clc
    close all
    a = [300:10:360 10:10:150];
    a = a + 5*(rand(size(a))-0.5);
    a(randperm(length(a),10))=-99;
else
    a = hd';
end
a(a == -99)= nan;
b = isnan(a);
c = [1 b]; c(end)=[];
d = [b 1]; d(1)=[];
e = (b + c + d)==3; % constriction
f = [false e]; f(end)=[];
g = [e false]; g(1)=[];
h = e | f | g; % dilation

% unwrap, fill missings and rewrap
k = unwrap(a/180*pi);
i = fillmissing(k,'spline');
l = mod(i,2*pi)*180/pi;

j= l.*~h+h*-99;
hd_filtered = j';
j(j==-99)=nan;
if nargout == 0
    plot(a,'+')
    hold on
    plot(j,'o')
end
end
