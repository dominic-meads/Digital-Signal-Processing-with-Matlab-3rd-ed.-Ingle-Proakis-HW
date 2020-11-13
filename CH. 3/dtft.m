function X = dtft(n,x,k)
% DESCRIPTION
%  Performs an approximate Discrete Time Fourier Transform 
%  the approximation is due to the fact that the output is 
%  not continous. ex. X = dtft(n,x,k); 
%  
% NOTE : Only for finite duration signals.
%
% INPUT VARIABLES
%  n = number of samples 
%  x = array containing time domain signal samples
%  k = evenly spaced divisions of omega (frequency)
%
% OUTPUT VARIABLES
%  X = the transformed output
%  Various graphs of signal aspects
%
% REFERENCES
%  Adapted from "Digital Signal Processing Using MATLAB 3rd ed." - Ingle, V.,
%  Proakis, J. (pg. 65). 
%
% DOCUMENTATION
%  ver 1.0 by Dominic Meads  5/8/2020
%  filename: dtft.m
%
% ENGINEER'S COMMENTS
%  The authors reccommend to use this more as an
%  excersize rather than a full function. As stated above, use the DFT for
%  better results.
%

kmax = length(k)-1;

w = (2*pi/kmax)*k;  % calculates the evenly spaced frequencies 
X = x * (exp(-j*2*pi/kmax)) .^(n'*k); % calculates the DTFT
end
