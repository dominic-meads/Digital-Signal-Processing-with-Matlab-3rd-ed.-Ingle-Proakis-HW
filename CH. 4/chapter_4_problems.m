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
