function [x,fval,exitflag] = optimization(K,M,N,B,F,scaleB,scaleF,p_max,R_k,g_k,N_0,C_k,E,E2,T2,options)

% Objective Function %
fun =@(x) obj(x);

% Initial Variable Guess %
numOfParameters = 3*K+1;
x0 = ones(numOfParameters,1);

% Variable Bounds %
lb = zeros(numOfParameters,1);
ub = [Inf; B*ones(M,1); F*ones(M,1); p_max*ones(M,1);...
    B*ones(N,1); F*ones(N,1); p_max*ones(N,1)];

% Linear Inequality Constraints %
[A,b] = LinearConstraints(M,N,B,F,scaleB,scaleF);

% Inequality Constraints %
Aeq = [];
beq = [];

% Scaling factor
scaling_factor = [1; scaleB*ones(M,1); scaleF*ones(M,1); ones(M,1);...
    scaleB*ones(N,1); scaleF*ones(N,1); ones(N,1)];

% Non-linear Inequality Constraints
nonlcon = NonLinearConstraints(R_k,g_k,N_0,C_k,E,E2,T2,M,N,K,...
    scaleB,scaleF);

[x_scaled,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
x = x_scaled.*scaling_factor;
end