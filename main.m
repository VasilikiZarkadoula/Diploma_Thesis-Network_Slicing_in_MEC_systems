% Convex Optimization Problem
close all;
clear;
clc;

% Constants %
K = 10;                             % total users
M = 5;                              % low latency users
N = K - M;                          % low energy consumption users

p_max = 0.5;                        % max transmission power of users
B = 10^6;                           % system total bandwidth
F = 10^11;                          % total computation capacity of cloud

E = 0.02:0.01:0.1;                  % energy threshold for low latency users
E2 = 0.007;                         % energy threshold for low energy users
T2 = 0.1;                           % latency threshold for low energy users

radius = 1000;                       % radius of a circular area
C_k = 1000*ones(K,1);               % CPU cycles
R_k = 1000*ones(K,1);               % bits
a = 2;                              % path-loss exponent
N_0 = 10^(-174 / 10);               % noise variance

L = 1;

fval = zeros(L,1);
exitflag = zeros(L,1);
count = zeros(1,length(E));
T_mean = zeros(1,length(E));
T_median = zeros(1,length(E));
for i = 1:length(E)
    for j = 1:L
        d_k = radius*rand(K,1);             % distance between user and BS
        h_k = exprnd(2,K,1);                % short-term fading
        g_k = (d_k.^(-a)).*(abs(h_k).^2);   % channel gain
        
        % Initial Variable Guess %
        B_m = 10^3 + (5*10^4 - 10^3)*rand(M,1);
        F_m = 10^8 + (5*10^9 - 10^8)*rand(M,1);
        p_m = 0.25*ones(M,1);
        
        B_n = 10^3 + (5*10^4 - 10^3)*rand(N,1);
        F_n = 10^8 + (5*10^9 - 10^8)*rand(N,1);
        p_n = 0.12*ones(N,1);
        
        T_0 = 0.1;
        
        
        x0 = [T_0; B_m; F_m; p_m; B_n; F_n; p_n];

        
        % Objective Function %
        fun = @(x) x(1);
        
        % Variable Bounds %
        lb = zeros(size(x0));  % all variables are positive
        ub = [Inf; B*ones(M,1); F*ones(M,1); p_max*ones(M,1);...
                B*ones(N,1); F*ones(N,1); p_max*ones(N,1)];
            


        % Linear Equality - Inequality Constraints %
        [A,b] = LinearConstraints(x0,M,N,K,p_max,B,F);
        
        Aeq = [];
        beq = [];
        
        % Non-linear Inequality Constraints 
        nonlcon = NonLinearConstraints(R_k,g_k,N_0,C_k,E(i),E2,T2,M,N);

        
        % Find optimal solution %
        options = optimoptions(@fmincon,'Display','iter',...
         'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-7,'StepTolerance',1e-16, ...
         'MaxFunctionEvaluations',1e+6,'MaxIterations', 1e+5);
        [x,fval(i),exitflag(i)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

        
%         if exitflag(j) == 0 || exitflag(j) == 2
%             fval(j) = 0;
%             count(i) = count(i) + 1;
%         end
       
    end 
    

end

% h1 = plot(E2,T_mean);
% hold on
% plot(E2,T_mean2,'r');
% xlabel('Energy threshold for low-latency users')
% ylabel('Latency threshold for low-latency users')
% title('p_max -> 0.05:0.1:10 | 500 iter');
% saveas(h1,'C:\Users\User\Desktop\Διπλωματική\Matlab\graphs\change_pmax\median\iter1000_5','jpeg');
% 
% figure;
% h2 = plot(p_max,T_mean);
% xlabel('Maximim transmission power of users')
% ylabel('Latency threshold for low-latency users')
% title('p_max -> 0.05:0.1:10 | 500 iter');
% saveas(h2,'C:\Users\User\Desktop\Διπλωματική\Matlab\graphs\change_pmax\mean\iter1000_5','jpeg');

% figure;
% plot()



