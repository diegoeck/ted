function [BkP] = experiment_Fu(M,W)
%This function computes the basis of the covariance matrix
%
%Consider the spectrum:
%
% Phi(omega) = sum_(i=1)^N c_i {delta(omega-omega_i) + delta(omega+omega_i)}
%
% Then
%
%inv(P) = sum_(i=1)^N c_i BkP{i} + R0
%
%Usage: [BkP] = experiment_Fu(M,W);
%
%               M: Experiment Design Model
%
%               W: Vector of Frequencies omega_i
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

for a=1:length(W)
    w=W(a);

    F=[];

    C1=freqresp(1+M.F.'*M.theta,w);
    C2=freqresp((M.B.'*M.theta)/(1+M.F.'*(M.theta))^2,w);
    C3=freqresp((1+M.D.'*M.theta)/(1+M.C.'*(M.theta)),w);
    
    for i=1:length(M.theta);
        F=[F;freqresp(M.B(i),w)/C1*C3-freqresp(M.F(i),w)*C2*C3];
    end

    BkP{a}=(F*F'+conj(F*F'))/(2*pi*M.sigma);
    
end


