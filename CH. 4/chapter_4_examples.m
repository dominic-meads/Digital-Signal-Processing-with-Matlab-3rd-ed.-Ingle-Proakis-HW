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
