function x = DblDiode(x,k)
%     x1 = x.^3;
%     x = x1./max(abs(x));
%% Model of Diode
% I = Isat exp(q/(kT) *V), 
% where k Boltsman's const, T temp in kelvin
% for T = 25 Celsius, q/(kT) = 40
if ~exist('k')
    k = 25e-3;%exp(-1);
end
x = x.*(1-Exp(-abs(x/k), 4)); %approximate to 4 terms, (
% x = x.*(1-exp(-abs(x/k)));
end%end function