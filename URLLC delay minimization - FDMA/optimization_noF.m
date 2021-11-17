function [x,fval,exitflag] = optimization_noF(K,U,M,B,scaleB,p_max,Fk,L,g,No,C,Eu,Em,Tm,options) 

% Objective Function %
fun =@(x) obj(x);

% Initial Variable Guess  %
numOfParameters = 2*K+1;
x0 = ones(numOfParameters,1);

% Variable Bounds %
lb = zeros(numOfParameters,1);
ub = [Inf;  B*ones(U,1); p_max*ones(U,1); B*ones(M,1); p_max*ones(M,1)];

% Linear Equality %
[A,b] = LinearConstraints_noF(U,M,B,scaleB);

% Inequality Constraints %
Aeq = [];
beq = [];

% Scaling factors
scaling_factor = [1; scaleB*ones(U,1); ones(U,1); scaleB*ones(M,1); ones(M,1)];

% Non-linear Inequality Constraints
nonlcon = NonLinearConstraints_noF(L,g,No,C,Fk,Eu,Em,Tm,U,M,K,scaleB);

[x_scaled,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
x =  x_scaled.*scaling_factor;
 end