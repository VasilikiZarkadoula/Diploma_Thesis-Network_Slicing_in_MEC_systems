% Convex Optimization Problem
close all;
clear;
clc;

% Constants %
K = 10;                             
M = 5;                              
N = K - M;                          

p_max = 0.5;                    
B = 1.8*10^6; 
F = 10^11;                         

E = 5*1e-3;                         
E2 = [1e-5 2.5*1e-5 5*1e-5 7.5*1e-5 1e-4 2.5*1e-4 5*1e-4 7.5*1e-4 1e-3];                       
T2 = 1;                           % latency threshold for low energy users

radius = 500;                       
C_k = 1000*ones(K,1);               
R_k = 500*ones(K,1);                
a = 2;                             
N_0 = 10^(-174 / 10);              

B_k = (B/K)*ones(K,1);
% F_k = (F/K)*ones(K,1);

L = 20;

fval = zeros(L,1);
exitflag = zeros(L,1);
count = zeros(1,length(E2));
T_mean = zeros(1,length(E2));

fval2 = zeros(L,1);
exitflag2 = zeros(L,1);
count2 = zeros(1,length(E2));
T_mean2 = zeros(1,length(E2));

for i = 1:length(E2)
    for j = 1:L
        d_k = radius*rand(K,1);
        h_k = exprnd(2,K,1);
        g_k = (d_k.^(-a)).*(abs(h_k).^2);
        
        % Initial Variable Guess %
        
        %         B_m = 10^3 + (5*10^4 - 10^3)*rand(M,1);
        B_m = [427313;320456;309044;241804;211595];
        F_m = 10^9 + (10^10 - 10^9)*rand(M,1);
        p_m = 0.2*ones(M,1);
%         p_m = p_max*rand(M,1);
        
        B_n = [2614;21834;24396;33152;207792];
        F_n = 10^7 + (5*10^8 - 10^7)*rand(N,1);
        p_n = 0.2*ones(N,1);
%         p_n = p_max*rand(N,1);

        F_m2 = 10^7 + (5*10^8 - 10^7)*rand(M,1);
        F_n2 = 10^9 + (10^10 - 10^9)*rand(M,1);
        
        T_0 = 0.01;
        
        x0 = [T_0; B_m; F_m; p_m; B_n; F_n; p_n];
        x0_2 = [1; F_m2; p_m; F_n2; p_n];
        %x0_3 = [T_0; B_m; p_m; B_n; p_n];
        
        % Objective Function %
        fun = @(x) x(1);
        
        % Variable Bounds %
        lb = zeros(size(x0));
        ub = [Inf; B*ones(M,1); F*ones(M,1); p_max*ones(M,1);...
            B*ones(N,1); F*ones(N,1); p_max*ones(N,1)];
        
        lb2 = zeros(size(x0_2));
        ub2 = [Inf; F*ones(M,1); p_max*ones(M,1); F*ones(N,1); p_max*ones(N,1)];
        
        % lb3 = zeros(size(x0_3));
        % ub3 = [Inf; B*ones(M,1); p_max*ones(M,1); B*ones(N,1); p_max*ones(N,1)];
        
        
        % Linear Equality - Inequality Constraints %
        [A,b] = LinearConstraints(M,N,B,F);
        [A2,b2] = LinearConstraints_noB(M,N,F);
        % [A3,b3] = LinearConstraints_noF(x0_3,M,N,K,p_max,B);
        
        Aeq = [];
        beq = [];
        
        % Non-linear Inequality Constraints
        nonlcon = NonLinearConstraints(R_k,g_k,N_0,C_k,E,E2(i),T2,M,N);
        nonlcon2 = NonLinearConstraints_noB(R_k,g_k,N_0,C_k,B_k,E,E2(i),T2,M,N);
        % nonlcon3 = NonLinearConstraints_noF(R_k,g_k,N_0,C_k,F_k,E(i),E2,T2,M,N);
        
        % Find optimal solution %
        options = optimoptions(@fmincon,...
            'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-6,'StepTolerance',1e-12, ...
            'MaxFunctionEvaluations',1e+6,'MaxIterations', 1e+5);
        [x,fval(j),exitflag(j)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        [x2,fval2(j),exitflag2(j)] = fmincon(fun,x0_2,A2,b2,Aeq,beq,lb2,ub2,nonlcon2,options);
        % [x3,fval3(j),exitflag3(j)] = fmincon(fun,x0_3,A3,b3,Aeq,beq,lb3,ub3,nonlcon3,options);
        
        if exitflag(j) == 0 || exitflag(j) == 2
            fval(j) = 0;
            count(i) = count(i) + 1;
        end
        
        if exitflag2(j) == 0 || exitflag2(j) == 2
            fval2(j) = 0;
            count2(i) = count2(i) + 1;
        end
        disp([i,j]);
    end
    newFval = fval(fval~=0);
    newFval2 = fval2(fval2~=0);
    T_mean(i) = mean(newFval);
    T_mean2(i) = mean(newFval2);
end


h1 = plot(E2,T_mean);
hold on
plot(E2,T_mean2,'r');
xlabel('Energy threshold for low-energy users')
ylabel('Latency threshold for low-latency users')
legend('variable bandwidth','equally shared bandwidth')
% saveas(h1,'C:\Users\User\Desktop\Διπλωματική\Matlab\graphs\change_bandwidth\iter100_E2_0.001_0.005_0.05','PNG');
% 
% figure;
% h2 = plot(p_max,T_mean);
% xlabel('Maximim transmission power of users')
% ylabel('Latency threshold for low-latency users')
% title('p_max -> 0.05:0.1:10 | 500 iter');
% saveas(h2,'C:\Users\User\Desktop\Διπλωματική\Matlab\graphs\change_pmax\mean\iter1000_5','jpeg');

% figure;
% plot()



