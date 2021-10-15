close all;
clear;
clc;

K = 16;       
S = K/2;        

p_max_u = (10^(24 / 10))/1000;      
p_max_m = (10^(21 / 10))/1000;  

BW = 10^6;
Bs = BW/S;            
L = 5e+4;
radius_m = 200;
radius_u = 500;   
dmin_k = 20;              
dmin_b = 30;             
a = 2;  
N0 = -174;                  
No = 10^((N0 - 30) / 10);       
                  
Tm = 1.5;
Eu = 5e-2;
Tu = 0.05;
b = L/(Eu*Bs);

iterations = 100;
count = 0;
e = 1e-20;

Pu = zeros(iterations,S);
Pm = zeros(iterations,S);
for i = 1:iterations
    
    flag = 0;   % break loop flag
    
    [dm,~,~] = distance_calc(dmin_k, dmin_b, radius_m , K/2);
    [du,~,~] = distance_calc(dmin_k, radius_m, radius_u , K/2);
    gu = matchURLLC(K,S,du,a,Bs,No);
    hm = randomGain(S);
    
    % define fixed P(z-1)
    pu_opt = p_max_u*ones(1,S);
    pm_opt = p_max_m*ones(1,S);
    
    do = true;
    z = 1;  % iteration index - while loop
    
    while do
        
        %%% step 1 - for fixed P(z-1) -> compute assignement %%%%
        pu_star = pu_opt;
        pm_star = pm_opt;
        % Em_opt_sum(z) = sum(pm_star.*Tm);
        
        [hung_matching, gm_hung] = matchMMTC_new(K,S,hm,gu,dm,pm_star,pu_star,L,a,Bs,No);
        if (length(gm_hung) ~= S)
            break;
        end
        
        %------------------------------------------------------------------
        
        %%%% step 2 - for obtained assignement -> compute P(z) %%%
        pu_opt = zeros(1,S);
        pm_opt = zeros(1,S);
        Em_opt = zeros(1,S);
        for s = 1:S
            aa = gu(s);
            b = L/(Eu*Bs);
            c = -(b* log(2)* 2^(b/aa))/(aa);
            root = (-b* log(2) - aa* lambertw(-1,c))/(aa* b* log(2));
            z1 = (1/gu(s))*(2^(L/(Bs*Tu)) -1);
            
            if z1 < real(root) && z1 < p_max_u
                pu_opt(s) = z1;
                pm_opt(s) = (1/gm_hung(s))*(2^(L/(Bs*Tm)) -1)*(pu_opt(s)*gu(s) + 1);
                % Em_opt(s) = pm_opt(s)*Tm;
            else
                % No feasible solution - break for loop
                count = count + 1;
                flag = 1;
                break;
            end
        end
        if(flag == 1) 
            % No feasible solution - break while loop
            break;
        end
        
        do = (any((abs(pu_opt - pu_star)) >= e)) || (any((abs(pm_opt - pm_star))  >= e));
        z = z+1;
        % Em_opt_sum(z) = sum(Em_opt);
    end
    Pu(i,:) = pu_opt;
    Pm(i,:) = pm_opt;
end








    
    

