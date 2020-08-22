function H = freqresp(b,a,w)
% DESCRIPTION
%  Calculates the frequency response of a LTI system described by a
%  difference equation.
%
% INPUT VARIABLES
%  b = numerator coefficient array (from standard form summation with x)
%  a = denominator coefficient array--a(1)=1--(from standard form summation with y)
%  w = frequency location array
%
% OUTPUT VARIABLES
%  H = Frequency response array
%
% REFERENCES
%  Adapted from "Digital Signal Processing Using MATLAB 3rd ed." - Ingle, V.,
%  Proakis, J. (pg. 101). 
%
% DOCUMENTATION
%  ver 1.0 by Dominic Meads  8/21/2020
%  filename: freqresp.m
%
% ENGINEER'S COMMENTS
%  see github (ch. 3 repo) for dtft function

a(1) = 1;   % change a(1) = 1 (see equation 3.21 on pg. 77)
num = 0;    % numerator
denom = 0;  % denominator

for m = 0:length(b)-1  % summation of numerator
    num = num + b(m+1)*exp(-j*m.*w);  % the 'b(m+1)' is because MATLAB cannot use zero as an index
end

for L = 0:length(a)-1  % summation of denominator
    denom = denom + a(L+1)*exp(-j*L.*w);  % the 'a(L+1)' is because MATLAB cannot use zero as an index
end

H = num./denom;  % final division (element by element) of arrays

end  % end function




