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
stem(n1,x1); axis([0 7 -4 4]);

n2 = 0:7;
x2 = [1 2 2 1 1 2 2 1]; 
%    ( x1[n])(x1[n-4])
figure;
stem(n2,x2); axis([0 7 -4 4]);

k = 0:200;

X1 = dtft(n1,x1,k);

X2 = dtft(n2,x2,k);

% check
w_check = (k*pi)/(length(k)-1);
X2_check = ones(1,201)*(X1 + exp(-4*j*w_check')*X1);
figure;
stem(abs(X2_check));
