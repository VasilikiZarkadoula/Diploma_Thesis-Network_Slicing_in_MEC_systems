function [x,fval,exitflag] = optimization_noB(K,U,M,F,scaleF,p_max,Bk,L,g,No,C,Eu,Em,Tm,options)

% Objective Function %
fun =@(x) obj(x);

% Initial Variable Guess  %
numOfParameters = 2*K+1;
x0 = ones(numOfParameters,1);

% Variable Bounds %
lb = zeros(numOfParameters,1);
ub = [Inf;  F*ones(U,1); p_max*ones(U,1); F*ones(M,1); p_max*ones(M,1)];

% Linear Inequality Constraints  %
[A,b] = LinearConstraints_noB(U,M,F,scaleF);

% Inequality Constraints %
Aeq = [];
beq = [];

% Scaling factors
scaling_factor = [1; scaleF*ones(U,1); ones(U,1);...
                     scaleF*ones(M,1); ones(M,1)];

% Non-linear Inequality Constraints
nonlcon = NonLinearConstraints_noB(L,g,No,C,Bk,Eu,Em,Tm,U,M,K,scaleF);

[x_scaled,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
x =  x_scaled.*scaling_factor;
end