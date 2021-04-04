% Convex Optimization Problem

clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % Constants %
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 30;                     % total users
M = 1:(K-1);                % low latency users
N = K - M;                  % low energy consumption users

rng default
R_k = 10^3*rand(K,1);       % bits
C_k = 1000*rand(K,1);       % CPU cycles
h = 10^(-3);                % channel gain
N_0 = 10^(-14);             % noise variance

p_max = 0.5;                % max transmission power of users
B = 10^6;                   % system total bandwidth
F = 10^11;                  % total computation capacity of cloud

E = 0.15;                   % energy threshold for low latency users
E2 = 0.01;                  % energy threshold for low energy users
T2 = 0.1;                   % latency threshold for low energy users



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Initial Variable Guess %
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fval = zeros(length(M),1);
for i=1:length(M)
    
    rng default
    B_m = 10^3 + (5*10^4 - 10^3)*rand(M(i),1);
    F_m = 10^8 + (5*10^9 - 10^8)*rand(M(i),1);
    p_m = repmat(0.25,[M(i),1]);
    
    B_n = 10^3 + (5*10^4 - 10^3)*rand(N(i),1);
    F_n = 10^8 + (5*10^9 - 10^8)*rand(N(i),1);
    p_n = repmat(0.12,[N(i),1]);
    
    T_0 = 0.09;
    
    x0 = [T_0; B_m; F_m; p_m; B_n; F_n; p_n];
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                         % Objective Function %
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fun = @(x) x(1);
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                        % Variable Bounds %
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lb = zeros(size(x0));  % all variables are positive
    ub = [];
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                % Linear Equality - Inequality Constraints %
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [A,b] = LinearConstraints(x0,M(i),N(i),K,p_max,B,F);
    
    Aeq = [];
    beq = [];
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                % Non-linear Inequality Constraints %
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nonlcon = NonLinearConstraints(R_k,h,N_0,C_k,E,E2,T2,M(i),N(i));
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                    % Find optimal solution %
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options = optimoptions(@fmincon,'StepTolerance',1e-18,...
                            'MaxFunEval',inf,'MaxIter',Inf);
    [x,fval(i),exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

for i = 1:length(M)
    fprintf('K = %i, M = %i, N = %i -> T = %0.3f ms\n',K,M(i),N(i),1000* fval(i));
end