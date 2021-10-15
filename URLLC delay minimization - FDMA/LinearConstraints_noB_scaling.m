function [A,b] = LinearConstraints_noB_scaling(M,N,F,scaleF)

A = [0, scaleF*ones(1,M), zeros(1,M), scaleF*ones(1,N), zeros(1,N)];
b = F;
