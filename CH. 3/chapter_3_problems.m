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
figure;
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
