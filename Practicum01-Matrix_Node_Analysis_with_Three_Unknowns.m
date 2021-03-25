clear all, close all, clc

G = 1/30*[10 -3; -3 9];

I_in = [6;9];
    
% GV=I;
V = inv(G)*I_in;

%% Output
disp('Conductance Matix is')
G

disp('Incoming Currents are')
I_in

disp('Node Voltages are')
V
