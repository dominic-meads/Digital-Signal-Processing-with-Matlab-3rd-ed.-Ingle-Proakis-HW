%% chapter 2 problems and notes
%% P2.1 
% generate the following sequences, plot using stem function
%
% * 1) 3*delta(n+2) + 2*delta(n) - delta(n-3) + 5*delta(n-7) -- -5 <= n <= 15
% * 2) 10*u(n) - 5*u(n-5) - 10*u(n-10) + 5*u(n-15) -- -5 <= n <= 15
% * 3) 5*cos(0.49*pi*n) + 5*cos(0.51*pi*n) -- -200 <= n <= 200

n1 = -5:15;
x1 = 3*impseq(-2,-5,15)+2*impseq(0,-5,15)-impseq(3,-5,15)+5*impseq(7,-5,15);
n2 = -5:15;
x2 = 10*stepseq(0,-5,15)-5*stepseq(5,-5,15)-10*stepseq(10,-5,15)+5*stepseq(15,-5,15);
n3 = -200:200;
x3 = 5*cos(0.49*pi*n3) + 5*cos(0.51*pi*n3);

figure;
subplot(2,2,1); stem(n1,x1); title('plot of function 1');
subplot(2,2,2); stem(n2,x2); title('plot of function 2');
subplot(2,2,3); stem(n3,x3); title('plot of function 3');

%% P2.5
% * 1) generate plot of exp(0,1,pi*n) for -100 <= n <= 100 using 'stem'
%      (plot both real and imaginary parts).
% * 2) generate plot of cos(0.1*n) for -20 <= n <= 20 using 'stem'

n1 = -100:100;
x1 = exp(0.1*pi*n1);
n2 = -20:20;
x2 = cos(0.1*n2);

figure;
subplot(1,2,1); stem(n1,real(x1)); title('plot of function 1 (real)');
subplot(1,2,2); stem(n1,imag(x1)); title('plot of function 1 (imaginary)');

% note: no imaginary component, therefore graph of imaginary is zero. 

figure; stem(n2,x2); title('plot of function 2');

%% P2.7
% * 1). modify the evenodd function to accept any arbitray sequence and
%     decomposes it into conjugate-symmetric and conjugate-antisymmetric 
%     parts. (see evenodd_conj.m)
% * 2). decompose the following sequence and plot:
%       x[n] = 10*exp[(-0.1 + 0.2*pi*j)*n]    0 <= n <= 10

n = 0:10;
x = 10*exp((-0.1 + 0.2*pi*j)*n);
[xe, xo, m] = evenodd_conj(x,n);

figure;
subplot(2,2,1); stem(n,real(x)); title('signal real part'); axis([-10,10,-10,10]);
subplot(2,2,2); stem(m,real(xe)); title('conj-symm real part'); axis([-10,10,-10,10]);
subplot(2,2,3); stem(m,real(xo)); title('conj-antisymm real part'); axis([-10,10,-10,10]);

figure;
subplot(2,2,1); stem(n,imag(x)); title('signal imag part'); axis([-10,10,-10,10]);
subplot(2,2,2); stem(m,imag(xe)); title('conj-symm imag part'); axis([-10,10,-10,10]);
subplot(2,2,3); stem(m,imag(xo)); title('conj-antisymm imag part'); axis([-10,10,-10,10]);

% the figures clearly show how xe and xo sum to the original x function

%% P2.8 
% The operation of signal dilation (or decimation or down-sampling) is defined by
%     y(n) = x(nM)
% in which the sequence x(n) is down-sampled by an integer factor M. For example, if
%     x(n) = {. . . ,−2, 4,3,−6, 5,−1, 8, . . .}
%                          ^
%  the down-sampled sequences by a factor 2 are given by
%     y(n) = {. . . ,−2,3, 5, 8, . . .}
%                       ^
% * 1. Develop a MATLAB function dnsample that has the form
%           function [y,m] = dnsample(x,n,M)
%      Downsample sequence x(n) by a factor M to obtain y(m)
%      to implement the above operation. Use the indexing mechanism of MATLAB with
%      careful attention to the origin of the time axis n = 0.
% * 2. Generate x(n) = sin(0.125πn), − 50 ≤ n ≤ 50. Decimate x(n) by a factor of 4 to
%      generate y(n). Plot both x(n) and y(n) using subplot and comment on the results.
% * 3. Repeat the above using x(n) = sin(0.5πn), − 50 ≤ n ≤ 50. Qualitatively discuss the
%      effect of down-sampling on signals.

