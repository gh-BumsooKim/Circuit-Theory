clear all, close all, clc

% Overdamped Case
R = 50;
L = 0.03;
C = 2e-04;
t = (0:0.001:1)*0.05;
vc1 = func_sereis_RLC_unit_step_vc(R, L, C, t);

% Critically Damped Case
L = 0.0425;
C = 1e-4;
vc2 = func_sereis_RLC_unit_step_vc(R, L, C, t);

% Underdamped Case
L = 0.1;
C = 5e-5;
vc3 = func_sereis_RLC_unit_step_vc(R, L, C, t);

plot(t, vc1, 'b', t, vc2, 'r', t, vc3, 'g');

function vc = func_sereis_RLC_unit_step_vc(R, L, C, t)

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
