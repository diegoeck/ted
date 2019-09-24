function [C,M,t]=experiment_e_trace(BkP,R0,energia);
%
%   TRACE OF TEH COST FUNCTION!!! I NEED TO CORRECT THE HELP
%
%This function computes optimal input spectrum
%
%Consider that:
%
% Phi(omega) = sum_(i=1)^N c_i {delta(omega-omega_i) + delta(omega+omega_i)}
% inv(P) = sum_(i=1)^N c_i BkP{i}
%
%The function solves the problem
%
% max t
%
% c_i>0
% sum(c_i)<energy
% inv(P)>t
%
%Usage: [C,M,topt]=experiment_e(Fu,R0,energy);
%
%               Fu: Basis of covariance matrix
%
%               energy: energy of the input signal
%
%               R0: matrix XXX
%
%               C: optimal spectrum [c_1 c_2 ... c_N]
%
%               M: optimal covariance matrix (inv(P))
%
%               topt: optimal t
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox

% Initialization
Mu=size(BkP,2);
n=size(BkP{1},1);

setlmis([]);

% Definition of the variables
[C,nn,sC]=lmivar(2,[1 Mu]);
for i=1:Mu
    r(i)=lmivar(3,sC(i));
end
%t=lmivar(1,[1 1]);
t=lmivar(1,[n 1]);


% Maximum Variance
l=newlmi;   
for a=1:Mu
lmiterm([-l,1,1,r(a)],1,BkP{a});
end
lmiterm([-l,1,1,0],R0);
lmiterm([-l,1,2,0],eye(n));
lmiterm([-l,2,2,t],1,1);





%lmiterm([l,1,1,t],eye(n),1);
%lmiterm([l,1,1,t],1,1);


% Positive Spectrum
for a=1:Mu
l=newlmi;   
lmiterm([l,1,1,r(a)],-1,1);
end

% Minimize Energy
l=newlmi;   
for a=1:Mu
lmiterm([l,1,1,r(a)],2,1/(2*pi));
end
lmiterm([l,1,1,0],-energia);



LMIS=getlmis;
%
%lmiinfo(LMIS);
%
%  Definition of the criterion trace(X)+trace(Y)
%
v=decnbr(LMIS);
c=zeros(v,1);
%for j=1:v
%    c(j)=-(defcx(LMIS,j,t));
%end
for j=1:v
c(j) = sum(diag(defcx(LMIS,j,t)));
end


% Solve the optimization problem
options=[0.00000001,1000,-1,100,1];
[copt,xopt]=mincx(LMIS,c,options);
%[copt,xopt]=mincx(LMIS,c)


% Create Solution
C=dec2mat(LMIS,xopt,C);
t=dec2mat(LMIS,xopt,t);
M=R0;
for a=1:Mu
        M=M+C(a)*BkP{a};
end

