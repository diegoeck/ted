function [u]=experiment_realization_cos(r,tsinal);
%This function computes optimal input signal
%
%Usage: [u]=experiment_realization_cos(C,length);
%
%               C: optimal spectrum [c_1 c_2 ... c_N]
%
%               length: length of the input signal
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

    N=length(r);
    sinal=randn(1,tsinal);
    sinal=sinal-mean(sinal);
    sinal=sinal/std(sinal);
    u=filter(firminphase([fliplr(r(2:N))/2 r(1) r(2:N)/2 ]),1,sinal); 
    u=u';




