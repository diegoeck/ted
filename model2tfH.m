function [H0]=model2tfH(M0);
%This function converts a Experiment Design Model M to the associated
%transfer function H.
%
%Usage: [H]=model2tfH(M);
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

H0=(1+M0.C.'*(M0.theta))/(1+M0.D.'*(M0.theta));
H0=simptf(H0);

