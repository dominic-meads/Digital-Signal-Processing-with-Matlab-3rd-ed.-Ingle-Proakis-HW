%% chapter 2 problems and notes
%% P2.1 
% generate the following sequences, plot using stem function
%
% * 1) 3*delta(n+2) + 2*delta(n) - delta(n-3) + 5*delta(n-7) -- -5 <= n <= 15
% * 2) 10*u(n) - 5*u(n-5) - 10*u(n-10) + 5*u(n-15) -- -5 <= n <= 15
% * 3) 5*cos(0.49*pi*n) + 5*cos(0.51*pi*n) -- -200 <= n <= 200
% * 4) e^0.1*n[u(n+20) - u(n-10)]

n1 = -5:15;
x1 = 3*impseq(-2,-5,15)+2*impseq(0,-5,15)-impseq(3,-5,15)+5*impseq(7,-5,15);
n2 = -5:15;
x2 = 10*stepseq(0,-5,15)-5*stepseq(5,-5,15)-10*stepseq(10,-5,15)+5*stepseq(15,-5,15);
n3 = -200:200;
x3 = 5*cos(0.49*pi*n3) + 5*cos(0.51*pi*n3);
n4 = -50:50;
x4 = exp(0.1*n4).*(stepseq(-20,-50,50)-stepseq(10,-50,50));

figure;
subplot(2,2,1); stem(n1,x1); title('plot of function 1');
subplot(2,2,2); stem(n2,x2); title('plot of function 2');
subplot(2,2,3); stem(n3,x3); title('plot of function 3');
subplot(2,2,4); stem(n4,x4); title('plot of function 4');

%% P2.3
% Generate the folowoing periodic sequnces and plot their samples using the
% stem function over the indicated amount of periods

% 1). [...,-2,-1,0,1,2,...]  plot 5 periods
%                ^ 

seq1 = [-2,-1,0,1,2];  % sequence
nseq1 = -2:2;
N1 = 5; % period in number of samples
P1 = 5; % amount of periods to plot
n1 = nseq1(1):P1*N1-(abs(nseq1(1))+1);
x1 = [seq1,zeros(1,(P1-1)*N1)];

for k = 1+N1:P1*N1
    x1(k) = x1(k-N1);
end

figure;
subplot(1,2,1); stem(n1,x1); title('5 periods of sequence 1');






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
% the figures clearly show how xe and xo sum to the original x function


% interestingly, the graphs show (and example 3.6) that a conjugate symmetric
% sequence will have even symmetry on the magnitude/real part, and odd 
% symmetry on the angle/imag part
magX = abs(xe);  % for graphs
angX = angle(xe);
realX = real(xe);  % divide into real and imaginary parts
imagX = imag(xe);

figure('Color', [1 1 1]);
subplot(2,2,1); plot(m,magX); grid off;
 title('Magnitude part');
subplot(2,2,2); plot(m,realX); grid off;
 title('Real part');
subplot(2,2,3); plot(m,angX); grid off;
 title('Angle part');
subplot(2,2,4); plot(m,imagX); grid off;
 title('Imaginary part');


%% P2.8 
% The operation of signal dilation (or decimation or down-sampling) is defined by
%     y(n) = x(nM)
% in which the sequence x(n) is down-sampled by an integer factor M. For example, if
%     x(n) = {. . . ,2,4,3,6,5,1,8, . . .}
%                        ^
%  the down-sampled sequences by a factor 2 are given by
%     y(n) = {. . . ,2,3,5,8, . . .}
%                      ^
% * 1. Develop a MATLAB function dnsample that has the form
%           function [y,m] = dnsample(x,n,M)
%      Downsample sequence x(n) by a factor M to obtain y(m)
%      to implement the above operation. Use the indexing mechanism of MATLAB with
%      careful attention to the origin of the time axis n = 0.
% * 2. Generate x(n) = sin(0.125*pi*n),  -50 <= n <= 50. Decimate x(n) by a factor of 4 to
%      generate y(n). Plot both x(n) and y(n) using subplot and comment on the results.
% * 3. Repeat the above using x(n) = sin(0.5*pi*n),  -50 <= n <= 50. Qualitatively discuss the
%      effect of down-sampling on signals.



%% P2.13
% for x[n] = cos(0.2*pi*n)+0.5*cos(0.6*pi*n) and echo component aplha(n-k),
% find:
% * 1). the 



