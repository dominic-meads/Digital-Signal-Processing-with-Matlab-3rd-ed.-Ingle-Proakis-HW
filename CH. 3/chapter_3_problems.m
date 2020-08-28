%% P3.1
% see dtft.m

% 4). x[n] = [4,3,2,1,1,2,3,4] plot dtft and comment on angle graph
%             ^

x = [4 3 2 1 1 2 3 4];
n = 0:7;
k = 0:100;
w = (pi/100)*k;  % calculates the evenly spaced frequencies 

X = dtft(n,x,k);

magX = abs(X);  % for graphs
angX = angle(X);
realX = real(X);  % divide into real and imaginary parts
imagX = imag(X);

figure('Color', [1 1 1]);
sgtitle('Discrete Time Fourier Transform');
subplot(2,2,1); plot(w*0.5/pi,magX); grid off; axis tight;
xlabel('frequency in units of pi'); title('Magnitude');
subplot(2,2,2); plot(w*0.5/pi,realX); grid off; axis tight;
xlabel('frequency in units of pi'); title('Real part');
subplot(2,2,3); plot(w*0.5/pi,angX); grid off; axis tight;
xlabel('frequency in units of pi'); title('Angle');
subplot(2,2,4); plot(w*0.5/pi,imagX); grid off; axis tight;
xlabel('frequency in units of pi'); title('Imaginary part');

% the angle graph is odd-symmetric


%% P3.2
% 1). After solving mathematically, the answer is X2(e^jw) = X1(e^jw) + e^-4jw * X1(e^jw)
%     with the first term coming from 0<=n<=3, and the second term coming
%     from the property of sample-shifting in the dtft and 4<=n<=7
% 2). 

n1 = 0:3;
x1 = [1 2 2 1];
figure;
subplot(2,2,1); stem(n1,x1); title('x1'); axis([0 7 -4 4]); xlabel('sample [n]'); ylabel('x[n]');

n2 = 0:7;
x2 = [1 2 2 1 1 2 2 1]; 
%    ( x1[n])(x1[n-4])
subplot(2,2,2); stem(n2,x2); title('x2'); axis([0 7 -4 4]); xlabel('sample [n]'); ylabel('x[n]');

k = 0:200;

X1 = dtft(n1,x1,k);

X2 = dtft(n2,x2,k);

% check
w_check = (pi/max(k))*k;
X2_check = X1 + exp(-4*j*w_check).*X1;  % use element by element multiplication (acts as a scalar)
subplot(2,2,3); stem(w_check/pi,abs(X2_check)); title('My DTFT (analyzed mathematically'); xlabel('Frequency as a fraction of PI'); ylabel('Magnitude');
subplot(2,2,4); stem(w_check/pi,abs(X2)); title('MATLAB computed DTFT (computed)'); xlabel('Frequency as a fraction of PI'); ylabel('Magnitude');
error = max(abs(X2 - X2_check))  % nominal error ~0%


%% P3.3
% 1). Find the DTFT of x[n] = 2(0.5)^n * u[n+2]

%     After solving analytically, X(e^jw) = 8e^(2jw) + 4e^(jw) + (2e^jw)/(e^jw - 0.5)

%     The shift in the unit step function causes the range for the
%     summation to be -2 <= n <= inf. A geometric series (3rd term in the
%     X(e^jw) equation above), is only defined a(r)^n when |r| < 1, or for
%     all n > 0. However, since our n includes negative samples, the
%     equation must account for those (when n < 0, |r| >= 1). Once all 
%     n < 0 terms are calculated, a geometric series for n < 0 can be used 
%     since |r| is now < 1.    

n = -2:200;                 % the sequence is zero for all n < -2 (due to shift)
u = stepseq(-2,-2,200);     % generate shifted unit step function
x = 2 * power(0.5,n) .* u;  % use element by element mult for each sample
k = 0:500;
w = k*pi/max(k);

X = dtft(n,x,k);

figure;
subplot(2,2,1); plot(w/pi,abs(X)); title('MATLAB Answer (Mag)'); xlabel('w as a fraction of pi'); ylabel('Magnitude'); axis([0 1 0 20]);
subplot(2,2,2); plot(w/pi,angle(X)*180/pi); title('MATLAB Answer (Ang)'); xlabel('w as a fraction of pi'); ylabel('Phase (degrees)'); axis([0 1 -200 200]);

% check answer 
X_check = 8*exp(2*j*w) + 4*exp(j*w) + (2*exp(j*w))./(exp(j*w)-0.5);  % my answer calculated analytically
subplot(2,2,3); plot(w/pi,abs(X_check)); title('My Answer (Mag)'); xlabel('w as a fraction of pi'); ylabel('Magnitude'); axis([0 1 0 20]);
subplot(2,2,4); plot(w/pi,angle(X_check)*180/pi); title('My Answer (Ang)'); xlabel('w as a fraction of pi'); ylabel('Phase (degrees)'); axis([0 1 -200 200]);
error = max(abs(X-X_check))  % nominal error ~= 0%;


%% P3.5

% 1). 
n1 = -2:2;
x1 = [2 1 3 1 2];  % answer found mathematically/analytically 
figure;
subplot(2,1,1); stem(n1,x1); title('My sequence #1'); axis([-10 10 -5 5]);

% 2).
n2 = -2:8;
x2 = [4 0 0 0 -3 1 -3 0 0 0 4];
subplot(2,1,2); stem(n2,x2); title('My sequence #2'); axis([-10 10 -5 5]);


%% P3.6 

% 1). Time domain equation is sin(pi*n/3)/(pi*n) -- like a sinc function
% 5).     { cos(pi*[n-10])/[n-10], n != 0
%         { 0, n = 0

n = -10:10;
k = -500:500;
w = k*pi/max(k);

x_c = [-0.1 0.111111 -0.125 0.142857 -0.166666 0.2 -0.25 0.333333 -0.5 1 0 -1 0.5 -0.333333 0.25 -0.2 0.166666 -0.142857 0.125 -0.111111 0.1];  % calculated with calculator
x = cos(pi*(n))./(n);  % my answer found analytically
x(21) = 0; 
error = max(abs(x-x_c))  % nominal error

% show both dtfts match up
X_c = dtft(n,x_c,k);
X = dtft(n,x,k);
figure;
plot(w/pi,abs(X_c));
figure;
plot(w/pi,abs(X));

% note I forgot to add in time shift of -10 samples (n-10) but the analytic
% and mathematical answers would still match up with one another



%% P3.8

n = -25:25;
k = -500:500;
w = k*pi/max(k);
x = exp(j*0.1*pi*n).*(stepseq(0,-25,25) - stepseq(20,-25,25));  % time domain sequence
% split time domain sequence into real and imaginary parts
xr = real(x);
xi = j*imag(x);

X = dtft(n,x,k);
[Xe Xo m] = evenodd_conj(X,k);  % decompose into conjugate symmetric and conjugate-antisymmetric parts (m is the same as k for this case)
% the evenodd_conj function was developed in chapter 2 :)

figure; 
subplot(2,2,1); plot(w/pi,abs(X)); title('DTFT'); xlabel('omega as a fraction of pi'); ylabel('Magnitude');
subplot(2,2,2); plot(w/pi,angle(X)); title('DTFT'); xlabel('omega as a fraction of pi'); ylabel('Phase (radians)');
subplot(2,2,3); plot(w/pi,Xe); title('DTFT (REAL)'); xlabel('omega as a fraction of pi'); ylabel('Magnitude');
subplot(2,2,4); plot(w/pi,Xo); title('DTFT (IMAGINARY)'); xlabel('omega as a fraction of pi'); ylabel('Magnitude');

Xe_check = dtft(n,xr,k);  % take the DTFT of real(x) to see if it is equal to the conjugate-symmetric part of X(e^jw)
error_even = max(abs(Xe_check-Xe))  % 0 error

Xo_check = dtft(n,xi,k);  % take the DTFT of imaginary(x) to see if it is equal to the conjugate-antisymmetric part of X(e^jw)
error_odd = max(abs(Xo_check-Xo))  % 0 error

% The above shows that the DTFT of both xr(real x) and xi(j*imaginary x) is
% equal to the conjugate-symmetric and conjugate-antisymmetric components 
% (respectivley) of the DTFT of x. This also shows the opposite, if you 
% have a time domain series x, and take the DTFT of that (X), 
% by splitting X into its conjugate-symmetric and conjugate-antisymmetric
% parts, and taking the IDTFT of those induvidual parts, your result will
% be the real parts of the time domain sequence, and j times the imaginary
% parts of the time domain sequence (respectivley)





%% P3.16
% see freqresp.m


%% P3.17

% 1).  LPF with poor stopband attenuation (fc ~= .038)
b1 = [.2 .2 .2 .2 .2];
a1 = zeros(1,length(b1));  % no 'a' coefficients
k = 0:500;
w = k*pi/max(k);
H1 = freqresp(b1,a1,w);
figure;
subplot(1,2,1); plot(w/pi,abs(H1)); sgtitle('problem 1'); title('Frequency Response (Magnitude)'); xlabel('frequency as a fraction of pi'); ylabel('magnitude');
subplot(1,2,2); plot(w/pi,angle(H1)*180/pi); title('Frequency Response (angle)'); xlabel('frequency as a fraction of pi'); ylabel('phase (degrees)');

% 2).  Bandpass ( 0.24 <= fpassband <= 0.43)
b2 = [1 0 -2];
a2 = [0 -0.95 0.9025];
H2 = freqresp(b2,a2,w);
figure;
subplot(1,2,1); plot(w/pi,abs(H2)/max(abs(H2))); sgtitle('problem 2'); title('Frequency Response (Magnitude)'); xlabel('frequency as a fraction of pi'); ylabel('magnitude');
subplot(1,2,2); plot(w/pi,angle(H2)*180/pi); title('Frequency Response (angle)'); xlabel('frequency as a fraction of pi'); ylabel('phase (degrees)');




