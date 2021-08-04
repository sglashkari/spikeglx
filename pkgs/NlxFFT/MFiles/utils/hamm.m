function W = hamm(N)
% Calc hamming window
a=0.53836;
W = a - (1-a)*cos((2*pi*(0:(N-1)))/(N-1));
W = W(:);