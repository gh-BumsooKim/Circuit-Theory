%% Step Response of Series Dynamic RLC Circuit with Zero Initials
% Euler Method
clear all, close all, clc
f0 = 60;
zeta = 0.5;
R = 50;
V0 = 10;
NT0 = 5;
NTs = 1000;

%%
w0 = 2*pi*f0;
alpha = zeta*w0; % zeta = alpha/w0
L = R/(2*alpha); % alpha = R/(2*L);
C = 1/(w0^2*L); % w0^2 = 1/(LC);
T0 = 1/f0;
Ts = T0/NTs; % Sample Time
tmax = NT0*T0;
t = -T0:Ts:tmax;
N = length(t);
disp('--- R, L, C, alpha, zeta, tmax ---')
format short e
[R; L; C; alpha; zeta; tmax]

%% FDM Analysis
a = 1 + R*C/Ts + L*C/Ts^2;
b = R*C/Ts + 2*L*C/Ts^2;
c = -L*C/Ts^2;

xvC = zeros(1,N);
vs = zeros(1,N);
xi = zeros(1,N);
xvL = zeros(1,N);

ii1 = find(t>0);
vs(ii1) = V0;

for ii=3:N
    xvC(ii) = (b/a)*xvC(ii-1) + (c/a)*xvC(ii-2) + (1/a)*vs(ii);
    xi(ii) = C*(xvC(ii) - xvC(ii-1))/Ts;
    xvL(ii) = L*(xi(ii) - xi(ii-1))/Ts;
end

xvR = xi*R;
imax = max(xi);
tt1 = 0:T0/5:tmax;
vvc1 = func_series_RLC_unit_step_vc(R,L,C,tt1);
vvc = V0*vvc1;

%% Data Plot
subplot(3,1,1), plot(t,xi), %axis([-Tc max(t) -0.5 1.5])
grid on
xlabel('time(s)'), ylabel('i(t)')
axis([min(t) max(t) [-0.5 1.5]*imax]), legend('i(t)')

subplot(3,1,2), plot(t,vs,t,xvC,t,xvL,t,xvR,t,xvC+xvL+xvR)
axis([min(t) max(t) -5 V0+5]), grid on
xlabel('time(s)'), ylabel('Voltage')
legend('v_s','v_C','v_L','v_R','v_C+v_L+v_R')

subplot(3,1,3), plot(t,xvC,tt1,vvc,'+')
axis([min(t) max(t) -5 V0+5])
title('Vadilation Check !!!')
xlabel('time(s)'), ylabel('v_c(t)')
legend('Numerical','Analytic')

function vc = func_series_RLC_unit_step_vc(R, L, C, t)

alpha = R/(2*L);
w0 = 1/sqrt(L*C);
zeta = alpha/w0;

s1 = -zeta * w0 + w0 * sqrt(zeta^2 -1);
s2 = -zeta * w0 - w0 * sqrt(zeta^2 -1); 

A1 = s2/(s1-s2);
A2 = - 1 - A1;

% Overdamped
if (alpha > w0) 
    vc = 1 + A1*exp(s1*t) + A2*exp(s2*t);
    
% Critically Damped
elseif (alpha == w0)
    vc = 1 + (A1 + A2*t).*exp(s1*t); 
    
% Underdamped
elseif (alpha < w0)
    vc = 1 - (cos(w0*t) + alpha/w0*sin(w0*t)).*exp(-alpha*t);
end
    
end