function [G0]=model2tfG(M0);
%This function converts a Experiment Design Model M to the associated
%transfer function G.
%
%Usage: [G]=model2tfG(M);
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

G0=M0.B.'*M0.theta/(1+M0.F.'*(M0.theta));
G0=simptf(G0);

