close all;
clear;
clc;

% Globals %

K = 10;                     % total users
M = 5;                      % low latency users
N = K - M;                  % low energy consumption users

p_max = 0.5;                % max transmission power of users                  
B = 10^6;                   % system total bandwidth
F = 10^9;                   % total computation capacity of cloud

C_k = 100;                 % CPU cycles               
R_k = 500;                  % bits  
radius = 500;               % radius of a circular area       
a = 2;                      % path-loss exponent 
N0 = -174;                  % dBm/Hz
N_0 = 10^((N0 - 30) / 10);  % noise variance

E = 1e-3;                   % energy threshold for low-latency users
E2 = 1e-5;                  % energy threshold for low-energy users
T2 = 0.1;                   % latency threshold for low-energy users

B_k = B/K;                  % equally shared bandwidth case
F_k = F/K;                  % equally shared CPU frequency case

% Scaling Factros
scaleB = 1e+4;
scaleF = 1e+7;
scaleP = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective Function %
fun =@(x) obj_scaling(x);

% Inequality Constraints %
Aeq = [];
beq = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Basic optimization case %

% Initial Variable Guess %
numOfParameters = 3*K+1;
x0 = ones(numOfParameters,1);
% Variable Bounds %
lb = zeros(numOfParameters,1);
ub = [Inf; B*ones(M,1); F*ones(M,1); p_max*ones(M,1);...
    B*ones(N,1); F*ones(N,1); p_max*ones(N,1)];
% Linear Inequality Constraints %
[A,b] = LinearConstraints_scaling(M,N,B,F,scaleB,scaleF);
% Scaling factors
scaling_factor = [1; scaleB*ones(M,1); scaleF*ones(M,1); ones(M,1);...
    scaleB*ones(N,1); scaleF*ones(N,1); ones(N,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 

                    % Constant Bandwidth case %

% Initial Variable Guess  %
numOfParameters_2 = 2*K+1;
x0_2 = ones(numOfParameters_2,1);
% Variable Bounds %
lb_2 = zeros(numOfParameters_2,1);
ub_2 = [Inf;  F*ones(M,1); p_max*ones(M,1); F*ones(N,1); p_max*ones(N,1)];
% Linear Inequality Constraints  %
[A2,b2] = LinearConstraints_noB_scaling(M,N,F,scaleF);
% Scaling factors
scaling_factor_2 = [1; scaleF*ones(M,1); ones(M,1);...
                       scaleF*ones(N,1); ones(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Constant Frequency case %
                   
% Initial Variable Guess %
numOfParameters_3 = 2*K+1;
x0_3 = ones(numOfParameters_3,1);
% Variable Bounds %
lb_3 = zeros(numOfParameters_3,1);
ub_3 = [Inf;  B*ones(M,1); p_max*ones(M,1); B*ones(N,1); p_max*ones(N,1)];
% Linear Equality %
[A3,b3] = LinearConstraints_noF_scaling(M,N,B,scaleB);
% Scaling factors
scaling_factor_3 = [1; scaleB*ones(M,1); ones(M,1);...
                       scaleB*ones(N,1); ones(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L = 1;
fval = zeros(1,length(E2));
exitflag = zeros(1,length(E2));
newFval = zeros(L,length(E2));
% x_scaled_2 = zeros(numOfParameters_2,length(E2));
% x2 = zeros(numOfParameters_2,length(E2));
fval_2 = zeros(1,length(E2));
exitflag_2 = zeros(1,length(E2));
newFval_2 = zeros(L,length(E2));
% x_scaled_3 = zeros(numOfParameters_2,length(E2));
% x3 = zeros(numOfParameters_2,length(E2));
% fval_3 = zeros(1,length(E2));
% exitflag_3 = zeros(1,length(E2));
for j = 1:L
    
    g_k = channelGain(radius,a,K);
    
    for i = 1:length(E2)
        
        % Find optimal solution %
        options = optimoptions(@fmincon,...
            'SpecifyObjectiveGradient',true,...
            'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-20,...
            'StepTolerance',1e-30, 'MaxFunctionEvaluations',1e+5,...
            'MaxIterations', 1e+5);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Basic optimization case %
                    
        % Non-linear Inequality Constraints
        nonlcon = NonLinearConstraints_scaling(R_k,g_k,N_0,C_k,E,E2(i),...
            T2,M,N,K,scaleB,scaleF,scaleP);

        [x_scaled(:,i),fval(i),exitflag(i)] = fmincon(fun,x0,A,b,Aeq,beq,lb,...
            ub,nonlcon,options);
        x(:,i) =  x_scaled(:,i).*scaling_factor;
%         if exitflag(i) ~= 1 
%             fval(i) = NaN;
%         end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%                      % Constant Bandwidth case %
%         
        % Non-linear Inequality Constraints 
        nonlcon_2 = NonLinearConstraints_noB_scaling(R_k,g_k,N_0,C_k,...
            B_k,E,E2(i),T2,M,N,K,scaleF,scaleP);
        
        [x_scaled_2(:,i),fval_2(i),exitflag_2(i)] = fmincon(fun,x0_2,A2,b2,...
            Aeq,beq,lb_2,ub_2,nonlcon_2,options);
        x2(:,i) =  x_scaled_2(:,i).*scaling_factor_2;
%         if exitflag_2(i) ~= 1 
%             fval_2(i) = NaN;
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                     % Constant Frequency case %
                            
        % Non-linear Inequality Constraints
        nonlcon_3 = NonLinearConstraints_noF_scaling(R_k,g_k,N_0,C_k,...
            F_k,E,E2(i),T2,M,N,K,scaleB,scaleP);
        
        [x_scaled_3,fval_3(i),exitflag_3(i)] = fmincon(fun,x0_3,A3,b3,...
            Aeq,beq,lb_3,ub_3,nonlcon_3,options);
        x3(:,i) =  x_scaled_3.*scaling_factor_3;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([j,i]);
    end
    
    newFval(j,:) = fval;
    newFval_2(j,:) = fval_2;
%     newFval_3(j,:) = fval_3;

end
% T_mean(i) = mean(newFval);
% TT = mean(newFval,'omitnan');
plot(E2,fval,'-*');
% hold on
% plot(E2,fval_2,'-*');
temp = newFval;
temp(any(isnan(temp), 2), :) = [];
T_mean = mean(temp);
% temp2 = newFval_2;
% temp2(any(isnan(temp2), 2), :) = [];
% T_mean2 = mean(temp2);
% h1 = plot(E2,1e+3.*T_mean,'-*');
% hold on
% plot(E2,1e+3.*T_mean2,'-*');
% xlabel('Energy threshold for low-energy users (J)')
% ylabel('Latency threshold for low-latency users (msec)')
% legend('variable bandwidth','equally shared bandwidth','Location','best');
% grid on
% saveas(h1,'C:\Users\User\Desktop\Διπλωματική\Matlab\graphs\change_LL_users_Energy\iter100','jpeg');
% hold on
% plot(E2,fval_2,'-*');
% hold on
% plot(E2,fval_3,'-*');
% legend('variable bandwidth and CPU frequency','equally shared bandwidth','equally shared CPU frequency','Location','best');

% sumBm = sum(x(2:6,:));
% sumBn = sum(x(17:21,:));
sumFm = sum(x(7:11,:));
sumFn = sum(x(22:26,:));

sumFm2 = sum(x2(2:6,:));
sumFn2 = sum(x2(12:16,:));
% disp(sumBm);
% disp(sumBn);
disp(sumFm);
disp(sumFn);

disp(sumFm2);
disp(sumFn2);