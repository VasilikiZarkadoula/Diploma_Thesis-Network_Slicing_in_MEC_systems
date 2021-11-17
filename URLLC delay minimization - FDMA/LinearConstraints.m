function [A,b] = LinearConstraints(U,M,B,F,scaleB,scaleF)

A1 = [0,scaleB*ones(1,U),zeros(1,U),zeros(1,U),...
        scaleB*ones(1,M),zeros(1,M),zeros(1,M)];
b1 = B;

A2 = [0,zeros(1,U),scaleF*ones(1,U),zeros(1,U),...
        zeros(1,M),scaleF*ones(1,M),zeros(1,M)];
b2 = F;

A = [A1; A2];
b = [b1; b2];
