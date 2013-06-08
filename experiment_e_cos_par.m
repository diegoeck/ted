function [rropt,M,topt]=experiment_e_cos_par(BkP,R0,energia);
%This function computes optimal input spectrum
%
%It uses the Partial Correlation Approach
%
%Consider that:
%
% Phi(omega) = sum_(i=1)^(N) c_i * cos( omega * (i-1) ) + EXPANSION TO
% INFINITY
% inv(P) = R0 + sum_(i=1)^(N) c_i BkP{i} + EXPANSION TO INFINITY
%
%The function solves the problem
%
% max t
%
% c_i>0
% sum(c_i)<energy
% inv(P)>t
%
%Usage: [C,M,topt]=experiment_e_cos_par(BkP,R0,energia);
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

for i=1:Mu
r(i) =lmivar(1,[1 1]);
end
t =lmivar(1,[1 1]);

R =lmivar(3,r);

%T=lmivar(3,toeplitz([r(1) 2*r(2:Mu)]))


for l=1:Mu
    for m=l:Mu

        if l==m
        lmiterm([1,l,m,r(1+m-l)],-1,2);
        else
        lmiterm([1,l,m,r(1+m-l)],-1,1);
        
        end
    end
end

%lmiterm([1,1,1,r(1)],-1,1);
%lmiterm([1,2,2,r(1)],-1,1);
%lmiterm([1,2,1,r(2)],-1,1);

%lmiterm([1,1,1,T],-1,1);
%lmiterm([1,1,1,0],-1*eps);


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
options=[1e-9,1000,-1,100,1];
[copt,xopt]=mincx(LMIS,c,options);
%[copt,xopt]=mincx(LMIS,c);
%[copt,xopt]=feasp(LMIS);


rropt=dec2mat(LMIS,xopt,R);

topt=dec2mat(LMIS,xopt,t);


M=R0;
for a=1:Mu
        M=M+rropt(a)*BkP{a};
end

