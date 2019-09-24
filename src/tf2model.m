function [M0]=tf2model(G0,H0,sigma);
%This function converts a model described by the transfer functions G(z)
%and H(z) to and Experiment Design Model M.
%
%Usage: [M]=tf2model(G,H,sigma)
%
%The model M has an Box-Jenkins structure:
%
%y(t)=G(z) u(t) + H(z) e(t)
%
%G(z,theta)=  B^T(z) theta      H(z,theta)= 1+ C^T(z) theta
%            ---------------                ---------------
%            1+ F^T(z) theta                1+ D^T(z) theta
%
%The vectors B(z), C(z), D(z) and F(z) are composed by polynomios.  
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox
[numG,denG]=tfdata(G0,'v');
[numH,denH]=tfdata(H0,'v');

[zG,pG,kG]=zpkdata(G0,'v');
delay=length(pG)-length(zG);

nG=length(numG);
nH=length(numH);

%M0.nB=nG-delay;
%M0.nF=nG-1;
%M0.nC=nH-1;
%M0.nD=nH-1;

iz=tf([0 1],[1 0],1);

M0.theta=[numG(1+delay:length(numG)) denG(2:length(denG)) numH(2:length(numH)) denH(2:length(denH))]';

M0.sigma=sigma;


M0.B=[];
for i=1+delay:nG
M0.B=[M0.B; iz^(i-1)];
end
for i=1:length(denG)+nH+nH-3
M0.B=[M0.B; tf(0,1,1)];
end

M0.F=[];
for i=1+delay:nG
M0.F=[M0.F; tf(0,1,1)];
end
for i=2:nG
M0.F=[M0.F; iz^(i-1)];
end
for i=1:nH+nH-2
M0.F=[M0.F; tf(0,1,1)];
end

M0.C=[];
for i=1+delay:nG+nG-1
M0.C=[M0.C; tf(0,1,1)];
end
for i=2:nH
M0.C=[M0.C; iz^(i-1)];
end
for i=1:nH-1
M0.C=[M0.C; tf(0,1,1)];
end

M0.D=[];
for i=1+delay:nG+nG-1+nH-1
M0.D=[M0.D; tf(0,1,1)];
end
for i=2:nH
M0.D=[M0.D; iz^(i-1)];
end




