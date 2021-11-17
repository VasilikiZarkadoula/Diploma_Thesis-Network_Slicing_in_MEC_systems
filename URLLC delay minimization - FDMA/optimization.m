function [x,fval,exitflag] = optimization(K,U,M,B,F,scaleB,scaleF,p_max,L,g,No,C_k,Eu,Em,Tm,options)

% Objective Function %
fun =@(x) obj(x);

% Initial Variable Guess %
numOfParameters = 3*K+1;
x0 = ones(numOfParameters,1);

% Variable Bounds %
lb = zeros(numOfParameters,1);
ub = [Inf; B*ones(U,1); F*ones(U,1); p_max*ones(U,1);...
           B*ones(M,1); F*ones(M,1); p_max*ones(M,1)];

% Linear Inequality Constraints %
[A,b] = LinearConstraints(U,M,B,F,scaleB,scaleF);

% Inequality Constraints %
Aeq = [];
beq = [];

% Scaling factor
scaling_factor = [1; scaleB*ones(U,1); scaleF*ones(U,1); ones(U,1);...
                     scaleB*ones(M,1); scaleF*ones(M,1); ones(M,1)];

% Non-linear Inequality Constraints
nonlcon = NonLinearConstraints(L,g,No,C_k,Eu,Em,Tm,U,M,K,scaleB,scaleF);

[x_scaled,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
x = x_scaled.*scaling_factor;
end