%% P4.1 
clear all;
clc;
% determine z transform and ROC of following expressions:
%
% 1).  x(n) = [3 2 1 -2 -3] (origin at x(n) = 1)
%
%       Finite-duration sequence: zero for n <n1 and  n > n2
%       For finite-duration sequences, ROC is entire z-plane
%       if n1 < 0, +inf is not in ROC, if n2 > 0, 0 is not in ROC
%       in this example, n1 = -2 -> -2 < 0: ROC doesnt include inf
%       also, n2 = 2 -> 2 > 0: ROC does include 0

%       X(z) = 3z^2 + 2z + 1 - 2z^-1 -3z^-2   ROC: 0 < |z| < +inf
%
% 2). x(n) = (0.8)^n * u(n-2) 
%
%       X(z) = (0.64z^-2)/(1-0.8z^-1)  ROC: 0.8 < |z| < inf
% 
%       MATLAB check: find the first 8 samples of x(n) corresonding to X(z)
%                     and compare to samples of actual x(n) from above
a2 = [1 -0.8 0];
b2 = [0 0 0.64];
[delta2,n2] = impseq(0,0,7);
x2 = filter(b2,a2,delta2);
[x2checkpart1,n2_1] = stepseq(2,0,7); 
x2check = (0.8.^n2) .* x2checkpart1;

error1 = max(abs(x2check-x2))  % no error, the z-transform above is the same

% 3). x(n) = (n+1)(3)^n * u(n) TODO ANSWER WRONG

a3 = [1 -6 9];
b3 = [3 0 0];
[delta3,n3] = impseq(0,0,7);
x3 = filter(b3,a3,delta3);
[x3checkpart1,n3_1] = stepseq(0,0,7); 
x3check = ((n3+1).*(3.^n3)).* x3checkpart1;
error1 = max(abs(x3check-x3))  % no error, the z-transform above is the same

%% P4.2
clear all;
clc;
% Given x(n) = (0.9)^n * cos(pi*n/4) * u(n)
% and
%               [  x(n/2), n = 0, +-2, +-4, ...
%       y(n) = [  0,      otherwise
%
% 1). show Y(z) = X(z^2)
%          n  = [-inf, ... -2,    -1,    0,    1,    2,  ... +inf]
%        x(n) = [-inf, ...x(-2), x(-1), x(0), x(1), x(2), ...+inf]
%        X(z) = [-inf, ...x(-2)*z^2, x(-1)*z^1, x(0), x(1)*z^-1, x(2)*z^-2, ...+inf]
% 
%        y(n) = [-inf, ...x(-1),  0  , x(0) ,  0  , x(1), ...+inf]
%        Y(z) = [-inf, ...x(-1)*z^2, x(0), x(1)*z^-2, ...+int] = X(z^2)
%
% 2). Determine Y(z)
%
%                             1 - 0.6346z^-2
%        Y(z) = X(z^2) = -------------------------   |z| > 0.9
%                        1 - 1.2728z^-2 + 0.81z^-4
%
% 3). Verify Y(z)
%      (i will only be using non-negative n because the sequence is causal)
% create x(n)
n = 0:100; 
x = ((0.9).^n).*cos(pi*n/4);

% create y(n)
y = zeros(1,101);  
for k = 1:101
    if (mod(k,2)==0)  % n is even
       y(k) = x(k/2);
    else 
        y(k) = 0;
    end
end
y = circshift(y,-1);  
% the reason I have to do the shift is to make y match up with the value of
% n, not the index number (which cannot be zero in MATLAB)

% verify calculated X(z), take inverse z-transform and check with x 
bx = [1 -0.6364 0];
ax = [1 -1.2728 0.81];
[deltax,nx] = impseq(0,0,100);
x_check = filter(bx,ax,deltax);
error_x = max(abs(x_check-x))

% verify calculated Y(z), take inverse z-transform and check with y 
by = [1 0 -0.6364 0 0];
ay = [1 0 -1.2728 0 0.81];
[deltay,ny] = impseq(0,0,100);
y_check = filter(by,ay,deltay);
error_y = max(abs(y_check-y))

%% P4.3
clear all;
clc;
%
% determine the z-transform of the following functions using properties of
% the z-transform and the table of z-transform pairs. Express X(z) as a
% rational function in z^-1. Verify results in MATLAB and plot poles and
% zeros

% 1). x(n) = 2*delta(n-2) + 3*u(n-3)  -- note: casaul by definition
%
%             2z^-2 - z^-3
%      X(z) = ------------  ROC |z| > 1
%               1 - z^-1

% verify z-transform by taking the inverse z-transform and checking the
% first 50 samples of x(n) from above, and the x(n) output from the inverse
% z-transform
b1 = [0 0 2 1];
a1 = [1 -1 0 0];
[delta1,n1] = impseq(0,0,49);
x1 = filter(b1,a1,delta1);
[x1_check1,ncheck11] = impseq(2,0,49);  % generate the delta function separatley
[x1_check2,ncheck12] = stepseq(3,0,49);  % generate the unit-step function separatley
x1_check = 2*x1_check1 + 3*x1_check2;
errorx1 = max(abs(x1_check-x1))
figure; 
subplot(2,2,1); zplane(b1,a1); title('P4.3.1 Poles and Zeros'); xlabel('Real'); ylabel('Imaginary');


% 2). x(n) = 3*(0.75)^n * cos(0.3*pi*n)*u(n) + 4*(0.75)^n * sin(0.3*pi*n)*u(n)
%
%                   3 + 1.0948z^-1
%      X(z) = ---------------------------    ROC |z| > 0.75
%             1 - 0.8817z^-1 + 0.5625z^-2

% verify z-transform by taking the inverse z-transform and checking the
% first 50 samples of x(n) from above, and the x(n) output from the inverse
% z-transform
b2 = [3 1.0948 0];
a2 = [1 -0.8817 0.5625];
[delta2,n2] = impseq(0,0,49);
x2 = filter(b2,a2,delta2);
[x2_check1,ncheck21] = stepseq(0,0,49);  % generate the unit-step function separatley
x2_check = 3*((0.75).^ncheck21).*cos(0.3*pi*ncheck21).*x2_check1 + 4*((0.75).^ncheck21).*sin(0.3*pi*ncheck21).*x2_check1;
errorx2 = max(abs(x2_check-x2))
% figure; stem(x2); hold on; stem(x2_check);  % error slightly noticable
subplot(2,2,2); zplane(b2,a2); title('P4.3.2 Poles and Zeros'); xlabel('Real'); ylabel('Imaginary');


% 3). I attempted it but failed :(
b3 = [0, 0.866, 0.0306, -1.5069, 2.8809, -2.8059, 0.81];
a3 = [1, -4.3641, 8.1177, -7.9641, 4.1177, -0.9, 0];
[delta3,n3] = impseq(0,0,49);
x3 = filter(b3,a3,delta3);
[x3_check1,ncheck31] = stepseq(0,0,49);  % generate the unit-step function separatley
[x3_check2,ncheck32] = stepseq(2,0,49);  % generate the unit-step function separatley
x3_check = ncheck31.*sin(pi*ncheck31./3).*x3_check1 + ((0.9).^ncheck31).*x3_check2;
errorx3 = max(abs(x3_check-x3))  %% wayyyyy to large 
subplot(2,2,3); zplane(b3,a3); title('P4.3.3 Poles and Zeros'); xlabel('Real'); ylabel('Imaginary');


%% P4.5
clear all;
clc;

% given Z{x[n]} = X(z) = 1/(1+0.5z^-1)   |z| >= 0.5
% find the z-transform of the following expressions

% 1). x1[n] = x[3-n] + x[n-3]

%               z^3 + 0.5z^2 + z^-3 + 0.5z^-4
%       X1(z) = -----------------------------     |z| >= 2
%                  0.5z + 1.25 + 0.5z^-1

% MATLAB check ERROR LINE 188 First denominator filter coefficient must be non-zero.
% first, find the original x[n] by numerically with the inverse z-transform
b = [1 0];
a = [1, 0.5];
[delta,n] = impseq(0,0,49);  
x = filter(b,a,delta);

% now make x1
x1 = zeros(1,50);  % assumes system is causal
for i = 1:50
    if (3-i <= 0)
        if (i-3 <= 0)
            x1 = 0;
        else
            x1(i) = 0 + x(i-3);
        end  % if (i-3...
    elseif (3-i > 0 && i-3 > 0)    
        x1(i) = x(3-i) + x(i-3);
    end  %if (3-i...
end  % for i = 1:50

% find inverse z-transform of X1(z) and test with x1[n]
b2 = [1 0.5 0 0 0 0 1 0.5];
a2 = [0 0 0.5 1.25 0.5 0 0 0];
[delta2,n2] = impseq(0,0,49);
x1_check = filter(b2,a2,delta2);
errorx1 = max(abs(x1_check-x1))
figure; zplane(b2,a2); title('P4.5.1 Poles and Zeros'); xlabel('Real'); ylabel('Imaginary');


%% P4.11.1
clear all;
clc;
%
% using partial fraction expansion, find the inverse z-transform of 
%
%                    1 - z^-1 -4z^-2 + 4z^-3
% 1). X1(z) = --------------------------------------   sequence right-sided
%             1 - (11/4)z^-1 + (13/8)z^-2 -(1/4)z^-3


b = [1 -1 -4 4];
a = [1 -11/4 13/8 -0.25];
[R p C] = residuez(b,a)

% x1[n] = u[n](-10(0.5)^n + 27(0.25)^n) - 16*delta[n]

% MATLAB check;
[delta,n] = impseq(0,0,49);
x1 = filter(b,a,delta);
[x1_check1,ncheck11] = impseq(0,0,49);  % generate the delta function separatley
[x1_check2,ncheck12] = stepseq(0,0,49);  % generate the unit-step function separatley
x1_check = x1_check2.*(-10*(0.5).^ncheck11 + 27*(0.25).^ncheck11) - 16.*x1_check1;
errorx1 = max(abs(x1_check-x1))  % ANSWER CORRECT


%% P4.15.1
clear all;
clc;
%
% Given h[n] = 5(1/4)^n * u[n]

% i). Find H(z)
%                      5
%         H(z) = ------------    |z| > 0.25
%                1 - 0.25z^-1

% ii). Find the difference equation representation
%
%         y[n] - 0.25*y[n-1] = 5*x[n]

% iii). Plot the poles and Zeros
b = [5 0];
a = [1 -0.25];
figure; zplane(b,a);

% iv). find the output y[n] if x[n] = (1/4)^n * u[n]
b1 = [5 0 0];
a1 = [1 -0.5 0.0625];
[R p C] = residuez(b1,a1)
%        y[n] = 5(1/4)^n * u[n]
%        the output is just scaled by 5

% MATLAB check for part iv
n = 0:7; 
x = (1/4).^n;
xic = 0;  % no intial conditions on x[n]
y = filter(b,a,x,xic)

ycheck = 5*(0.25).^n
















