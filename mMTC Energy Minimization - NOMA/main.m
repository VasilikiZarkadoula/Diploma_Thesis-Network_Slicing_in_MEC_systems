% Diploma Thesis - mMTC Sum Energy Minimization with NOMA
% Simulation Results
% Vasiliki Zarkadoula
close all;
clear;
clc;

K = 10;                                 % total users       
S = K/2;                                % subchannels       
p_max_u = (10^(24 / 10))/1000;          % max transmission power of URLLC      
p_max_m = (10^(21 / 10))/1000;          % max transmission power of mMTC   
BW = 10^6;                              % system total bandwidth
Bs = BW/S;                              % bandwidth of each subchannel            
L = 5e+4;                               % bits
radius = 500;                           % radius of a circular area       
d_min_users = 20;                       % minimum distance between users              
d_min_bs = 30;                          % minimum distance from BS              
a = 2;                                  % path-loss exponent   
N0 = -174;                              % AWGN noise (dBm/Hz)                  
No = 10^((N0 - 30) / 10);       
                  
Tm = 1.5;                               % latency threshold of mMTCs   
Tu = 0.05;                              % latency threshold of URLLCs 
Eu = 0.1;                               % energy threshold of URLLCs 



%%% Iterative Algorithm - Joint subchannel and power allocation %%%

iterations = 100;
count = 0;
e = 1e-15;

Pu = zeros(iterations,S);
Pm = zeros(iterations,S);
Em_hung = zeros(iterations,S);
Em_exh = zeros(iterations,S);
Em_rand = zeros(iterations,S);
for i = 1:iterations
    
    flag = 0;   % break loop flag
    
    % find distances of mMTCs and URLLCs from BS 
    [dm,~,~] = distance_calc(d_min_users, d_min_bs, radius , S);
    [du,~,~] = distance_calc(d_min_users, d_min_bs, radius , S);
    
    gu = matchURLLC(K,du,a,Bs,No); % matching of URLLCs to subchannels
    hm = randomGain(S); 
    
    % randomly matching mMTCs to subchannels
    [rand_matching_opt, gm_rand] = matchMMTC_rand(K,hm,dm,a,Bs,No);
    
    % define fixed Pu(0),Pm(0) and A(0)
    pu_opt = p_max_u*ones(1,S);
    pm_opt_hung = p_max_m*ones(1,S);
    hung_matching_opt = zeros(1,S);
    pm_opt_exh = p_max_m*ones(1,S);
    exh_matching_opt = zeros(1,S);
    
    do = true;
    z = 1;      % iteration index - while loop
    while do
        
        %%% STEP 1 : for fixed P(z-1) -> compute assignement %%%
       
        pu_star = pu_opt;
        
        pm_star_hung = pm_opt_hung;
        pm_star_exh = pm_opt_exh;
        
        hung_matching_star = hung_matching_opt;
        exh_matching_star = exh_matching_opt;
        
        % Flags defining whether to use the Hungarian algorithm or the exhaustive search algorithm
        hungMethod = 0;
        exhMethod = 1;
        
        % optimally matching mMTCs to subchannels
        [hung_matching_opt, gm_hung] = matchMMTC_opt(K,hm,gu,dm,pm_star_hung,pu_star,L,a,Bs,No,hungMethod);
        [exh_matching_opt, gm_exh] = matchMMTC_opt(K,hm,gu,dm,pm_star_exh,pu_star,L,a,Bs,No,exhMethod);
        
        % Error check
        if (length(gm_hung) ~= S) || (length(gm_exh) ~= S) || (length(gm_rand) ~= S)
            flag = 1;
            pu_opt = zeros(1,S);
            pm_opt_hung = zeros(1,S);
            pm_opt_exh = zeros(1,S);
            break;
        end
        
        %--------------------------------------------------------------
        
        %%% STEP 2 : for obtained assignement -> compute P(z) %%%
        
        pu_opt = zeros(1,S);
        pm_opt_hung = zeros(1,S);
        Em_opt_hung = zeros(1,S);
        pm_opt_exh = zeros(1,S);
        Em_opt_exh = zeros(1,S);
        pm_opt_rand = zeros(1,S);
        Em_opt_rand = zeros(1,S);
        for s = 1:S
            aa = gu(s);
            b = L/(Eu*Bs);
            c = -(b* log(2)* 2^(b/aa))/(aa);
            root = (-b* log(2) - aa* lambertw(-1,c))/(aa* b* log(2));
            z1 = (1/gu(s))*(2^(L/(Bs*Tu)) -1);
            
            if z1 < real(root) && z1 < p_max_u
                pu_opt(s) = z1;
                
                % Closed form solutions when the Hungarian, Exhaustive and
                % Random method is followed respectively
                pm_opt_hung(s) = (1/gm_hung(s))*(2^(L/(Bs*Tm)) -1)*(pu_opt(s)*gu(s) + 1);
                Em_opt_hung(s) = pm_opt_hung(s)*Tm;
                
                pm_opt_exh(s) = (1/gm_exh(s))*(2^(L/(Bs*Tm)) -1)*(pu_opt(s)*gu(s) + 1);
                Em_opt_exh(s) = pm_opt_exh(s)*Tm;
                
                pm_opt_rand(s) = (1/gm_rand(s))*(2^(L/(Bs*Tm)) -1)*(pu_opt(s)*gu(s) + 1);
                Em_opt_rand(s) = pm_opt_rand(s)*Tm;
                
                % Feasibility check
                if pm_opt_hung(s) > p_max_m || pm_opt_exh(s) > p_max_m || pm_opt_rand(s) > p_max_m
                    % No feasible solution - break for loop s = 1:S
                    count = count + 1;
                    flag = 1;
                    break;
                end
            else
                % No feasible solution - break for loop s = 1:S
                count = count + 1;
                flag = 1;
                break;
            end
        end
        if(flag == 1)
            % No feasible solution - break while loop
            break;
        end
        do = (any((abs(pu_opt - pu_star)) >= e)) || ...
            (any((abs(pm_opt_hung - pm_star_hung))  >= e)) || ...
            ~isequal(hung_matching_star,hung_matching_opt) || ...
            (any((abs(pm_opt_exh - pm_star_exh))  >= e)) || ...
            ~isequal(exh_matching_star,exh_matching_opt);
        z = z+1;
        if (z >= 100)
            break;
        end
        %Em_opt_sum(z-1) = sum(Em_opt_hung);
    end
    %         if(flag == 1)
    %             % No feasible solution - break for loop l = 1:length(L)
    %             break;
    %         end
    Pu(i,:) = pu_opt;
    Pm(i,:) = pm_opt_hung;
    Em_hung(i,:) = Em_opt_hung;
    Em_exh(i,:) = Em_opt_exh;
    Em_rand(i,:) = Em_opt_rand;
%     Em_sum_opt_hung = sum(Em_opt_hung);
%     Em_sum_opt_exh = sum(Em_opt_exh);
%     Em_sum_opt_rand = sum(Em_opt_rand);
%     Em_hung(i,l) = Em_sum_opt_hung;
%     Em_exh(i,l) = Em_sum_opt_exh;
%     Em_rand(i,l) = Em_sum_opt_rand;
    
end

% Delete rows with infeasible solutions
Pu(Pu == 0) = NaN;
Pu(any(isnan(Pu),2),:) = [];
Pm(Pm == 0) = NaN;
Pm(any(isnan(Pm),2),:) = [];
Em_hung(Em_hung == 0) = NaN;
Em_hung(any(isnan(Em_hung),2),:) = [];
Em_exh(Em_exh == 0) = NaN;
Em_exh(any(isnan(Em_exh),2),:) = [];
Em_rand(Em_rand == 0) = NaN;
Em_rand(any(isnan(Em_rand),2),:) = [];

% Bar Plot
Pu_mean = mean(Pu);
Pm_mean = mean(Pm);
y = [Pu_mean; Pm_mean];
bar(1:S,y)
xlabel('Subchannel index');
ylabel('Average mMTC/URLLC transmission power (Watt)')
title('Power allocation for S = 5 subchannels')
legend('URLLC','mMTC','Location','northwest')
set(gca, 'YScale', 'log');

