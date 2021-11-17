% Diploma Thesis - URLLC Delay Minimization with FDMA
% Simulation Results
% Vasiliki Zarkadoula
close all;
clear;
clc;

K = 10;                     % total users
U = K/2;                    % URLLC users
M = K - U;                  % mMTC users
p_max = 0.5;                % max transmission power of users                  
B = 10^6;                   % system total bandwidth
F = 10^10;                  % total computation capacity of edge cloud
C = 100;                    % CPU cycles    
L = 5e+4;                   % bits
radius = 500;               % radius of a circular area       
a = 2;                      % path-loss exponent                             
N0 = -174;                  % AWGN noise (dBm/Hz)
No = 10^((N0 - 30) / 10);  


Eu = 1e-4;                  % energy threshold of URLLCs        		
Em = 1e-6;                  % energy threshold of mMTCs            
Tm = 0.3:0.1:1;             % latency threshold of mMTCs   

Bk = B/K;                   % equally shared bandwidth case
Fk = F/K;                   % equally shared CPU frequency case

% Scaling Factros
scaleB = 1e+4;
scaleF = 1e+8;

% Options for fmincon
options = optimoptions(@fmincon,'Display','off',...
    'SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-8,...
    'ConstraintTolerance',1e-20,'StepTolerance',1e-30,...
    'MaxFunctionEvaluations',1e+5,'MaxIterations', 1e+5);

% Find optimal solution - Monte Carlo Simulations %
iterations = 1000;

Tu_tot = zeros(iterations,length(Tm));
Tu_tot_2 = zeros(iterations,length(Tm));
Tu_tot_3 = zeros(iterations,length(Tm));
for j = 1:iterations
    g = channelGain(radius,a,K);
    for i = 1:length(Tm)
        
        % proposed optimization - dynamic bandwidth and CPU allocation
        [x,Tu,exitflag] = optimization(K,U,M,B,F,scaleB,scaleF,p_max,L,g,No,C,Eu,Em,Tm(i),options);
        % fixed bandwidth allocation case
        [x2,Tu_2,exitflag2] = optimization_noB(K,U,M,F,scaleF,p_max,Bk,L,g,No,C,Eu,Em,Tm(i),options);
        % fixed CPU capacity allocation case
        [x3,Tu_3,exitflag3] = optimization_noF(K,U,M,B,scaleB,p_max,Fk,L,g,No,C,Eu,Em,Tm(i),options);
        
        % check for feasibility of optimization
        if exitflag ~= 1
            Tu = NaN;
        end
        if exitflag2 ~= 1
            Tu_2 = NaN;
        end
        if exitflag3 ~= 1
            Tu_3 = NaN;
        end
        
        Tu_tot(j,i) = Tu;
        Tu_tot_2(j,i) = Tu_2;
        Tu_tot_3(j,i) = Tu_3;
    end
    disp(j);
end

% Remove iterations with non feasible solutions
Tu_tot(any(isnan(Tu_tot),2),:) = [];
Tu_tot_2(any(isnan(Tu_tot_2),2),:) = [];
Tu_tot_3(any(isnan(Tu_tot_3),2),:) = [];

% Display results
Tu_mean = mean(Tu_tot);
Tu_2_mean = mean(Tu_tot_2);
Tu_3_mean = mean(Tu_tot_3);
plot(Tm,Tu_mean,'-*');
hold on
plot(Tm,Tu_2_mean,'-*');
hold on
plot(Tm,Tu_3_mean,'-*');
xlabel('Latency threshold for mMTC users (sec)');
ylabel('Average Latency threshold for URLLC users (sec)')
legend('Proposed Resource Allocation','Fixed Bandwidth','Fixed CPU Frequency','Location','best')
hold off

perc2 = ((Tu_2_mean-Tu_mean)/Tu_mean)*100;
perc3 = ((Tu_3_mean-Tu_mean)/Tu_mean)*100;