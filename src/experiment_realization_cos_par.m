function [u]=experiment_realization_cos_par(r,tsinal);
%This function computes optimal input signal
%
%Usage: [u]=experiment_realization_cos_par(C,length);
%
%               C: optimal spectrum [c_1 c_2 ... c_N]
%
%               length: length of the input signal
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

N = length(r);
sinal = randn(1,tsinal);

r0 = r(1);
r1n = r(2:N)/2;

% Youla Walker Equation
A = toeplitz([ r0 r1n(1:(N-2)) ]);
f = A\(r1n');
k1 = r0-f'*r1n';

% Check stability of AR model
%abs(roots([1;-f]'))

u = filter([(k1)^(1/2)],[1;-f]',sinal); 
u = u';

    


