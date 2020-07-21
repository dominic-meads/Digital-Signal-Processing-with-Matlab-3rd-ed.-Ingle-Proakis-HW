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

w = (pi/100)*k;  % calculates the evenly spaced frequencies 
X = x * (exp(-j*pi/100)) .^(n'*k); % calculates the DTFT

magX = abs(X);  % for graphs
angX = angle(X);
realX = real(X);  % divide into real and imaginary parts
imagX = imag(X);

figure('Color', [1 1 1]);
subplot(2,2,1); plot(w/pi,magX); grid off;
xlabel('frequency in units of pi'); title('Magnitude part');
subplot(2,2,2); plot(w/pi,realX); grid off;
xlabel('frequency in units of pi'); title('Real part');
subplot(2,2,3); plot(w/pi,angX); grid off;function X = dtft(n,x,k)
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

w = (pi/100)*k;  % calculates the evenly spaced frequencies 
X = x * (exp(-j*pi/100)) .^(n'*k); % calculates the DTFT

end

xlabel('frequency in units of pi'); title('Angle part');
subplot(2,2,4); plot(w/pi,imagX); grid off;
xlabel('frequency in units of pi'); title('Imaginary part');
