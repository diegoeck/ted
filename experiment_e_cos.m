function [rropt,M,topt]=experiment_e_cos(BkP,R0,energia);
%This function computes optimal input spectrum
%
%Consider that:
%
% Phi(omega) = sum_(i=1)^(N) c_i * cos( omega * (i-1) )
% inv(P) = sum_(i=1)^(N) c_i BkP{i} + R0
%
%The function solves the problem
%
% max t
%
% c_i>0
% sum(c_i)<energy
% inv(P)>t
%
%Usage: [C,M,topt]=experiment_e_cos(Fu,R0,energy);
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

Mu=size(BkP,2);

n=size(BkP{1});
setlmis([]);

% Definition of the variables
t =lmivar(1,[1 1]);
Q =lmivar(1,[Mu-1 1]);
A = [zeros(1,Mu-2) 0; eye(Mu-2) zeros(Mu-2,1)];
B = [1 zeros(1,Mu-2)]';
[C,nn,sC]=lmivar(2,[1 Mu-1]);
[D,nn,sD]=lmivar(2,[1 1]);

r(1)=lmivar(3,sD(1));
for i=1:Mu-1
    r(i+1)=lmivar(3,sC(i));
end



lmiterm([-1,1,1,Q],1,1);
lmiterm([-1,1,1,Q],-A',A);
lmiterm([-1,1,2,Q],-A',B);
lmiterm([-1,2,2,Q],-B',B);
 
lmiterm([-1,1,1,0],zeros(Mu-1,Mu-1));
lmiterm([-1,1,2,-C],1,1);
lmiterm([-1,2,2,D],1,1);
lmiterm([-1,2,2,-D],1,1);


% Maximum Variance
l=newlmi;   
for a=1:Mu
lmiterm([-l,1,1,r(a)],1,BkP{a});
end
lmiterm([-l,1,1,0],R0);
lmiterm([l,1,1,t],eye(n),1);


% Minimize Energy
l=newlmi;   
lmiterm([l,1,1,r(1)],1,1);
lmiterm([l,1,1,0],-energia);



LMIS=getlmis;
%
%lmiinfo(LMIS);
%
%  Definition of the criterion trace(X)+trace(Y)
%
v=decnbr(LMIS);
c=zeros(v,1);
for j=1:v
    c(j)=-defcx(LMIS,j,t);
end


%
% Solve the optimization problem
%
options=[0.00000001,1000,-1,100,1];
[copt,xopt]=mincx(LMIS,c,options);
%[copt,xopt]=mincx(LMIS,c);
%[copt,xopt]=feasp(LMIS);


Qopt=dec2mat(LMIS,xopt,Q);
Copt=dec2mat(LMIS,xopt,C);
Dopt=dec2mat(LMIS,xopt,D);

topt=dec2mat(LMIS,xopt,t);

rropt = [Dopt Copt];

M=R0;
for a=1:Mu
        M=M+rropt(a)*BkP{a};
end
