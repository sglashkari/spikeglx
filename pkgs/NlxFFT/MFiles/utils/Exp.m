function y = Exp(x, N)
% y = Exp(x, N)
% approx exp(x) with N terms
% for error < 1%
% x = 1 2 3  4  5  6  7  8  9  10
% N = 5 7 9 10 12 13 15 16 18  19
% N ~ round(1.54x + 3.93)
if nargin <2
    N = 5;
end
    
n = 0:ceil(N-1);
for i = 1:numel(x)
    y(i) = sum((x(i).^n)./factorial(n));
end
y = reshape(y, size(x));
