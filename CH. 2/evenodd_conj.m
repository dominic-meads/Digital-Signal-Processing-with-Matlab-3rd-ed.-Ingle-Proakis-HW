function [xe, xo, m] = evenodd_conj(x,n)

% Real signal decomposition into conjugate-symmetric and conjugate-antisymmetric parts
% conjugate-symmetric -> xe(n) = xe*(-n)
% conjugate-antisymmetri -> xo(n) = -xo*(-n)
% 
% adapted from evenodd, edited by Dominic Meads
% -------------------------------------------------

% [xe, xo, m] = evenodd_conj(x,n)

%

m = -fliplr(n);

m1 = min([m,n]); m2 = max([m,n]); m = m1:m2;

nm = n(1)-m(1); n1 = 1:length(n);

x1 = zeros(1,length(m));

x1(n1+nm) = x; x = x1;

xe = 0.5*(x + conj(fliplr(x)));

xo = 0.5*(x - conj(fliplr(x)));
