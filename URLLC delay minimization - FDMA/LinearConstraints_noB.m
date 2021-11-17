function [A,b] = LinearConstraints_noB(U,M,F,scaleF)

A = [0, scaleF*ones(1,U), zeros(1,U), scaleF*ones(1,M), zeros(1,M)];
b = F;
