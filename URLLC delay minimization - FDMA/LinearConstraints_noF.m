function [A,b] = LinearConstraints_noF(U,M,B,scaleB)

A = [0, scaleB*ones(1,U), zeros(1,U), scaleB*ones(1,M), zeros(1,M)];
b = B;
