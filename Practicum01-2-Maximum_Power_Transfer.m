clear all, close all, clc

Vs = 10;
Rs = 50;
RL = 0:1:200;

%%
VL = Vs*RL./(Rs+RL);
IL = Vs./(Rs+RL);
PL1 = Vs^2*RL/(RL+Rs).^2;
PL = VL.*IL;

figure('name','Maximum Power Transfer')
subplot(3,1,1), plot(RL, VL)
xlabel('RL (Ohm)'), ylabel('VL (V)')
subplot(3,1,2), plot(RL, IL)
subplot(3,1,3), plot(RL, PL)

PLmax = max(PL);
iimax = find(PL == PLmax);
RLopt = RL(iimax)
