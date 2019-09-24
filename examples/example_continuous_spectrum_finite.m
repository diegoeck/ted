clc
clear all

G=(zpk([0.7],[0.8 0.9],1,1))
H=tf([1],[1],1);
sigma=1;

M0=tf2model(G,H,sigma);

energy=1;
length=1000;

N=30;
[BkP] = experiment_Fu_cos_ss(M0,N);
[R0] = experiment_R0_cos_ss(M0);

[r,invP,topt]=experiment_e_cos(BkP,R0,energy)

[u]=experiment_realization_cos(r,length);

W=0:0.001:3.14;
sinalrr=0;
for a=1:30
sinalrr=sinalrr+r(a)*cos((a-1)*W);
end
semilogx(W,sinalrr, 'b');

figure(2)
plot(u)

