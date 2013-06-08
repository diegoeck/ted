clc
clear all

G=(zpk([0.7],[0.8 0.9],1,1))
H=tf([1],[1],1);
sigma=1;

M0=tf2model(G,H,sigma);

energy=1;
length=1000;

Ini=3;
Freq=30;
W1=[];
for a=1:Freq
    W1(a)=10^(((a-1)*Ini/Freq)-Ini)*pi;
end

[BkP] = experiment_Fu(M0,W1);
[R0] = experiment_R0_cos_ss(M0);

[r,invP,topt]=experiment_e(BkP,R0,energy)

u=experiment_realization(r,W1,length);

figure(1)
semilogx( W1,r,'bo');

figure(2)
plot(u)

