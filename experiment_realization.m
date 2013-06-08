function [u]=experiment_realization(r,W,tsinal);
%This function computes optimal input signal
%
%Usage: [u]=experiment_realization(C,W,length);
%
%               C: optimal spectrum [c_1 c_2 ... c_N]
%
%               W: frequencies
%
%               length: length of the input signal
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

u=zeros(1,tsinal);
quantos=0;

for a=1:length(W)
%        if r(a)>max(r)/100
            quantos=quantos+1;
            u=u+sqrt(2*r(a)/pi)*sin(rand*2*pi+W(a)*(0:(tsinal-1)));
%        end
end
u=u';

