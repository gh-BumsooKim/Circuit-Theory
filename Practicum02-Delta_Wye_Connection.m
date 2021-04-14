clear all, close all, clc

syms Ra Rb Rc Rab Rbc Rca

% ABCD for T-network
%
%----Ra-------Rc---
%         |
%         Rb
%         |
%------------------

TT = [1 Ra; 0 1] * [1 0; 1/Rb 1] * [1 Rc; 0 1];

% ABCD for pi-network
%
%--------Rab--------
%    |        |
%    Rca      Rbc
%    |        |
%------------------

PP = [1 0; 1/Rca 1] * [1 Rab; 0 1] * [1 0; 1/Rbc 1];

xRb = 1/PP(2,1)
xRc = (PP(2,2) - 1) * xRb
xRa = (PP(1,1) - 1) * xRb

xRab = TT(1,2)
xRbc = xRab/(TT(1,1) - 1)
xRca = xRab/(TT(2,2) - 1)