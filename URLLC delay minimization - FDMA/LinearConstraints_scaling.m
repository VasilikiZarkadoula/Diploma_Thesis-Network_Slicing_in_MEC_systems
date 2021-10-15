function [A,b] = LinearConstraints_scaling(M,N,B,F,scaleB,scaleF)

A1 = [0,scaleB*ones(1,M),zeros(1,M),zeros(1,M),scaleB*ones(1,N),zeros(1,N),zeros(1,N)];
b1 = B;

A2 = [0,zeros(1,M),scaleF*ones(1,M),zeros(1,M),zeros(1,N),scaleF*ones(1,N),zeros(1,N)];
b2 = F;

A = [A1; A2];
b = [b1; b2];
