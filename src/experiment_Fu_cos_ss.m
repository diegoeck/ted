function [BkP] = experiment_Fu_cos_ss(M,N)
%This function computes the basis of the covariance matrix.
%It uses state-space representation for fast computetion.
%
%Consider the spectrum:
%
% Phi(omega) = sum_(i=1)^(N) c_i cos( omega * (i-1) )
%
% Then
%
%inv(P) = sum_(i=1)^(N) c_i BkP{i} + R0
%
%Usage: [BkP] = experiment_Fu(M,N);
%
%               M: Experiment Design Model
%
%               N: Size of the basis
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox
w=sym('w','real');

Fu=(1+M.D.'*(M.theta))/(1+M.C.'*(M.theta))*(M.B/(1+M.F.'*(M.theta))-(M.B.'*M.theta)*(M.F)/(1+M.F.'*(M.theta))^2);

Fuss=ss(Fu);

P=dlyap(Fuss.a,Fuss.b*Fuss.b');

BkP{1}=(Fuss.c*P*Fuss.c'+Fuss.d*Fuss.d')/(M.sigma);

for a=2:N
   
    Mp=Fuss.c*(Fuss.a)^(a-2)*(Fuss.a*P*Fuss.c'+Fuss.b*Fuss.d');
    BkP{a}=(Mp+Mp')/(2*M.sigma);
end


