%% ex. 4.9
% Compute the inverse z-transform of the following function
%                    
%  X(z) =  ___________1_____________   |z| > 0.9
%          (1-0.9z^-1)^2 (1+0.9z^-1)

% 1). find polynomial coeficients

b = 1;
a = poly([0.9,0.9,-0.9]);  % poly returns the coefficients of the polynomial given the roots

% 2). Find poles, residues, and direct terms (direct terms are a result of
%     partial fraction decomp when the order of the numerator exceeds or is
%     equal to the order of the denominator

[R,p,C] = residuez(b,a)

% now we have 
%                X(z) = ___0.25___ + ___0.25____ + ______0.5____
%                       (1+0.9z^-1)  (1-0.9z^-1)   (1-0.9z^-1)^2
%
% two repeated poles of 0.9 result in first and second order denominators

% ANSWER: analytically, the inverse z-transform is
%         x[n] = 0.25(-0.9)^n*u[n] + 0.75(0.9)^n*u[n] + 0.5*n*(0.9)^n*u[n]


%% Ex. 4.11
figure;

% b and a must be same length
b = [1 0]; 
a = [1 -0.9];

zplane(b,a);  % note input arguments for zplane are ROW VECTORS

[H w] = freqz(b,a,100,'whole'); % in matlab, the second half (bottom part) of unit circle starts at pi, so without "whole," w only runs to 0.99*pi
magH = abs(H); 
phaH = angle(H);
subplot(2,1,1); plot(w/pi,magH); grid on; xlabel('frequency in pi units'); ylabel('Magnitude'); title('Magnitude Response'); axis([0 1 0 10]);
subplot(2,1,2); plot(w/pi,phaH/pi); grid on; xlabel('frequency in pi units'); ylabel('Phase in pi units'); title('Phase Response'); axis([0 1 -.4 .4]);

%% Ex. 4.13
figure;

b = [1 0 -1];
a = [1 0 -0.81];
[R p C] = residuez(b,a)

k = 0:500;
w = k*pi/500;
H = freqz(b,a,w);
plot(w/pi,abs(H));
figure; plot(w/pi,angle(H));

%% Ex. 4.15

% solve the difference equation:
% 
%   y(n) = 1/3[x(n) + x(n-1) + x(n-2)] + 0.95*y(n-1) - 0.9025*y(n-2)  n>=0
% 
% with the initial conditions
% 
%   y(-1) = -2
%   y(-2) = -3
%   x(-1) = 1
%   x(-2) = 1
% and
%
%   x(n) = cos(pi*n/3)*u(n)
%

% time domain coefficients
b = [1 1 1]/3;
a = [1 -0.95 0.9025];

% initial conditions
X = [1 1];
Y= [-2 -3];
xic = filtic(b,a,Y,X)  % filtic is included in signal processing toolbox and calculates xic(z)

% X(z) z transform coefficients (one-sided z transform)
bxplus = [1 -0.5];
axplus = [1 -1 1];

% convolution can be used to multiply two NON-CAUSAL z-domain polynomials
% get Y(z) coefficients
ayplus = conv(a,axplus)
byplus = conv(b,bxplus) + conv(xic,axplus) % add the conv(xic,axplus) becuase of initial condition term

% calulate function using Partial Fraction Expansion
[R,p,C] = residuez(byplus,ayplus)

% find respective magnitudes and phases
Mp = abs(p)
Ap = angle(p)/pi

