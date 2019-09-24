clc
clear all

G = zpk([],[.8 .9],1,1);
H = zpk([0.75],[0.85],1,1);
sigma = 0.01;

M0 = tf2model(G,H,sigma);

energy = 1.5;
length = 1000;

N = 9;
[BkP] = experiment_Fu_cos_ss(M0,N);
[R0] = experiment_R0_cos_ss(M0);

[r,invP,topt] = experiment_e_cos_par(BkP,R0,energy);

[u] = experiment_realization_cos_par(r,length);

figure(1)
W=0:0.001:3.14;
sinalrr=0;
for a=1:N
sinalrr=sinalrr+r(a)*cos((a-1)*W);
end
semilogx(W,sinalrr, 'b');

figure(2)
plot(u)

figure(3)
plot(W,sinalrr)

figure(4)
pwelch(u)