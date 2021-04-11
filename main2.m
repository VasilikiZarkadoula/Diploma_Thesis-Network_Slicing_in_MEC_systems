% Convex Optimization Problem
close all;
clear;
clc;

% Constants %
K = 10;                             % total users
M = 5;                              % low latency users
N = K - M;                          % low energy consumption users

p_max = 0.2;                        % max transmission power of users
B = 2.5*10^5;                           % system total bandwidth
F = 10^11;                          % total computation capacity of cloud

E = 0.1:0.1:0.9;                  % energy threshold for low latency users
% E2 = 0.007;                         % energy threshold for low energy users
E2 = 0.001;
T2 = 0.5;                           % latency threshold for low energy users

radius = 500;                       % radius of a circular area
C_k = 1000*ones(K,1);               % CPU cycles
R_k = 500*ones(K,1);               % bits
a = 2;                              % path-loss exponent
N_0 = 10^(-174 / 10);               % noise variance

B_k = (B/K)*ones(K,1);

L = 100;

fval = zeros(L,1);
exitflag = zeros(L,1);
count = zeros(1,length(E));
T_mean = zeros(1,length(E));
fval2 = zeros(L,1);
exitflag2 = zeros(L,1);
T_mean2 = zeros(1,length(E));
for i = 1:length(E)
    for j = 1:L
        d_k = radius*rand(K,1);             % distance between user and BS
        h_k = exprnd(0.5,K,1);                % short-term fading
        g_k = (d_k.^(-a)).*(abs(h_k).^2);   % channel gain
        
        % Initial Variable Guess %
        B_m = 10^3 + (5*10^4 - 10^3)*rand(M,1);
        F_m = 10^8 + (5*10^9 - 10^8)*rand(M,1);
%         p_m = 0.25*ones(M,1);
        p_m = p_max*rand(M,1);
        
        B_n = 10^3 + (5*10^4 - 10^3)*rand(N,1);
        F_n = 10^8 + (5*10^9 - 10^8)*rand(N,1);
%         p_n = 0.12*ones(N,1);
        p_n = p_max*rand(N,1);
        
        T_0 = 0.1;
        
        
        x0 = [T_0; B_m; F_m; p_m; B_n; F_n; p_n];
        x0_2 = [T_0; F_m; p_m; F_n; p_n];
        
        % Objective Function %
        fun = @(x) x(1);
        
        % Variable Bounds %
        lb = zeros(size(x0));  % all variables are positive
        ub = [Inf; B*ones(M,1); F*ones(M,1); p_max*ones(M,1);...
                B*ones(N,1); F*ones(N,1); p_max*ones(N,1)];
            
        lb2 = zeros(size(x0_2));  
        ub2 = [Inf; F*ones(M,1); p_max*ones(M,1); F*ones(N,1); p_max*ones(N,1)];    


        % Linear Equality - Inequality Constraints %
        [A,b] = LinearConstraints(x0,M,N,K,p_max,B,F);
        [A2,b2] = LinearConstraints_noB(x0_2,M,N,K,p_max,F);
        
        Aeq = [];
        beq = [];
        
        % Non-linear Inequality Constraints 
        nonlcon = NonLinearConstraints(R_k,g_k,N_0,C_k,E(i),E2,T2,M,N);
        nonlcon2 = NonLinearConstraints_noB(R_k,g_k,N_0,C_k,B_k,E(i),E2,T2,M,N);
        
        % Find optimal solution %
        options = optimoptions(@fmincon,...
         'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-9,'StepTolerance',1e-14, ...
         'MaxFunctionEvaluations',1e+6,'MaxIterations', 1e+5);
        [x,fval(j),exitflag(j)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        [x2,fval2(j),exitflag2(j)] = fmincon(fun,x0_2,A2,b2,Aeq,beq,lb2,ub2,nonlcon2,options);
        
        if exitflag(j) == 0 || exitflag(j) == 2
            fval(j) = 0;
            count(i) = count(i) + 1;
        end
       
    end 
    newFval = fval(fval~=0);
    T_mean(i) = mean(fval);
    T_mean2(i) = mean(fval2);
end

h1 = plot(E,T_mean);
hold on
plot(E,T_mean2,'r');
xlabel('Energy threshold for low-energy users')
ylabel('Latency threshold for low-latency users')
legend('variable bandwidth','equally shared bandwidth')
% saveas(h1,'C:\Users\User\Desktop\ÄéðëùìáôéêÞ\Matlab\graphs\change_pmax\median\iter1000_5','jpeg');
% 
% figure;
% h2 = plot(p_max,T_mean);
% xlabel('Maximim transmission power of users')
% ylabel('Latency threshold for low-latency users')
% title('p_max -> 0.05:0.1:10 | 500 iter');
% saveas(h2,'C:\Users\User\Desktop\ÄéðëùìáôéêÞ\Matlab\graphs\change_pmax\mean\iter1000_5','jpeg');

% figure;
% plot()



