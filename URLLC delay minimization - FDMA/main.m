close all;
clear;
clc;

K = 10;                 % total users
M = K/2;                    % low latency users
N = K - M;                  % low energy consumption users

p_max = 0.5;                % max transmission power of users                  
B = 10^6;                   % system total bandwidth
F = 10^10;                  % total computation capacity of cloud

C_k = 100;                  % CPU cycles    
R_k = 5e+4;
radius = 500;               % radius of a circular area       
a = 2;                      % path-loss exponent                             
N0 = -174;                  % dBm/Hz
N_0 = 10^((N0 - 30) / 10);  % noise variance


E = 1e-4;         		
E2 = 1e-7:1e-7:1e-6;            
T2 = 0.1; 

B_k = B/K;                  % equally shared bandwidth case
F_k = F/K;                  % equally shared CPU frequency case

% Scaling Factros
scaleB = 1e+4;
scaleF = 1e+8;

% Options for fmincon
options = optimoptions(@fmincon,...
    'Display','off',...
    'SpecifyObjectiveGradient',true,...
    'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-20,...
    'StepTolerance',1e-30, 'MaxFunctionEvaluations',1e+5,...
    'MaxIterations', 1e+5);

% Find optimal solution %
L = 1000;

fval = zeros(1,length(E2));
exitflag = zeros(1,length(E2));
newFval = zeros(L,length(E2));
fval2 = zeros(1,length(E2));
exitflag2 = zeros(1,length(E2));
newFval2 = zeros(L,length(E2));
fval3 = zeros(1,length(E2));
exitflag3 = zeros(1,length(E2));
newFval3 = zeros(L,length(E2));
for j = 1:L
    
    g_k = channelGain(radius,a,K);
    
    for i = 1:length(E2)
        
        [x,fval(i),exitflag(i)] = optimization(K,M,N,B,F,scaleB,scaleF,p_max,R_k,g_k,N_0,C_k,E,E2(i),T2,options);
        [x2,fval2(i),exitflag2(i)] = optimization_noB(K,M,N,F,scaleF,p_max,B_k,R_k,g_k,N_0,C_k,E,E2(i),T2,options);
        [x3,fval3(i),exitflag3(i)] = optimization_noF(K,M,N,B,scaleB,p_max,F_k,R_k,g_k,N_0,C_k,E,E2(i),T2,options);

       
        if exitflag(i) ~= 1 
            fval(i) = NaN;
            
        end
        if exitflag2(i) ~= 1 
            fval2(i) = NaN;
            
        end
        if exitflag3(i) ~= 1 
            fval3(i) = NaN;
            
        end

        
    end
    disp(j);
    newFval(j,:) = fval;
    newFval2(j,:) = fval2;
    newFval3(j,:) = fval3;
end
% plot(1e-3*R_k,fval); hold on; plot(1e-3*R_k,fval2); hold on; plot(1e-3*R_k,fval3)

newFval(any(isnan(newFval),2),:) = [];
newFval2(any(isnan(newFval2),2),:) = [];
newFval3(any(isnan(newFval3),2),:) = [];

T = mean(newFval,'omitnan');
T_2 = mean(newFval2,'omitnan');
T_3 = mean(newFval3,'omitnan');
plot(1e-6*E2,1e+3*T,'-*');
hold on
plot(1e-6*E2,1e+3*T_2,'-*');
hold on
plot(1e-6*E2,1e+3*T_3,'-*');
xlabel('Energy threshold for mMTC users (μJ)');
ylabel('Average Latency threshold for URLLC users (msec)')
legend('Proposed Resource Allocation','Fixed Bandwidth','Fixed CPU Frequency','Location','best')
hold off

perc2 = ((T_2-T)/T)*100;
perc3 = ((T_3-T)/T)*100;


% sumBm = sum(x(2:6,:));
% sumBn = sum(x(17:21,:));
% disp(sumBm);
% disp(sumBn);
% sumFm = sum(x(7:11,:));
% sumFn = sum(x(22:26,:));
% disp(sumFm);
% disp(sumFn);
% sumFm2 = sum(x2(2:6,:));
% sumFn2 = sum(x2(12:16,:));
% disp(sumFm2);
% disp(sumFn2);


