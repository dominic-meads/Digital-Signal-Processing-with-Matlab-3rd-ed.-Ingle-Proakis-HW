%% P4.1
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
x2 = filter(b2,a2,delta2)
[x2checkpart1,n2_1] = stepseq(2,0,7); 
x2check = (0.8.^n2) .* x2checkpart1

error1 = max(abs(x2check-x2))  % no error, the z-transform above is the same

% 3). x(n) = (n+1)(3)^n * u(n) TODO ANSWER WRONG

a3 = [1 -6 9];
b3 = [3 0 0];
[delta3,n3] = impseq(0,0,7);
x3 = filter(b3,a3,delta3)
[x3checkpart1,n3_1] = stepseq(0,0,7); 
x3check = ((n3+1).*(3.^n3)).* x3checkpart1
error1 = max(abs(x3check-x3))  % no error, the z-transform above is the same

%% P4.2
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

% % verify calculated Y(z), take inverse z-transform and check with y 
by = [1 0 -0.6364 0 0];
ay = [1 0 -1.2728 0 0.81];
[deltay,ny] = impseq(0,0,100);
y_check = filter(by,ay,deltay);
error_y = max(abs(y_check-y))
