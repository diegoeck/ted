function [R0] = experiment_R0_cos_ss(M)
%This function computes matrix R0 used to compute invP.
%It uses state-space representation for fast computetion.
%
%Consider the spectrum:
%
% Phi(omega) = sum_(i=1)^(N) c_i * cos( omega * (i-1) )
%
% Then
%
%inv(P) = sum_(i=1)^(N) c_i BkP{i} + R0
%
%Usage: [R0] = experiment_R0_cos_ss(M);
%
%               M: Experiment Design Model
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox
w=sym('w','real');

Fu=(1+M.D.'*(M.theta))/(1+M.C.'*(M.theta))*(M.C/(1+M.D.'*(M.theta))-(1+M.C.'*M.theta)*(M.D)/(1+M.D.'*(M.theta))^2);

Fuss=ss(Fu);

P=dlyap(Fuss.a,Fuss.b*Fuss.b');

R0=(Fuss.c*P*Fuss.c'+Fuss.d*Fuss.d');



