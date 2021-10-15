function [x,fval,exitflag] = optimization_noB(K,M,N,F,scaleF,p_max,B_k,R_k,g_k,N_0,C_k,E,E2,T2,options)

% Objective Function %
fun =@(x) obj(x);

% Initial Variable Guess  %
numOfParameters = 2*K+1;
x0 = ones(numOfParameters,1);

% Variable Bounds %
lb = zeros(numOfParameters,1);
ub = [Inf;  F*ones(M,1); p_max*ones(M,1); F*ones(N,1); p_max*ones(N,1)];

% Linear Inequality Constraints  %
[A,b] = LinearConstraints_noB(M,N,F,scaleF);

% Inequality Constraints %
Aeq = [];
beq = [];

% Scaling factors
scaling_factor = [1; scaleF*ones(M,1); ones(M,1);...
    scaleF*ones(N,1); ones(N,1)];

% Non-linear Inequality Constraints
nonlcon = NonLinearConstraints_noB(R_k,g_k,N_0,C_k,...
    B_k,E,E2,T2,M,N,K,scaleF);

[x_scaled,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
x =  x_scaled.*scaling_factor;
end