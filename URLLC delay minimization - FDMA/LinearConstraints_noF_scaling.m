function [A,b] = LinearConstraints_noF_scaling(M,N,B,scaleB)

A = [0, scaleB*ones(1,M), zeros(1,M), scaleB*ones(1,N), zeros(1,N)];
b = B;
