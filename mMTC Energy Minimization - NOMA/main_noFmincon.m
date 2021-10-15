close all;
clear;
clc;

K = 8:2:30;       
S = K/2;        

p_max_u = (10^(24 / 10))/1000;      
p_max_m = (10^(21 / 10))/1000;  

BW = 10^6;
% Bs = BW/S;            
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
Tu = 0.03;

% b = L/(Eu*Bs);

iterations = 5000;

tic;

count = 0;


Em_sum_hung = zeros(iterations,length(K));
Em_sum_rand = zeros(iterations,length(K));
for k = 1:length(K)
    
    Bs = BW/S(k);
    for i = 1:iterations
        
        [dm,~,~] = distance_calc(dmin_k, dmin_b, radius_m , K(k)/2);
        [du,~,~] = distance_calc(dmin_k, radius_m, radius_u , K(k)/2);
        
        gu = matchURLLC(K(k),S(k),du,a,Bs,No);
        matchMMTC(K(k),S(k),gu,dm,p_max_m,p_max_u,L,a,Bs,No);
        [gm_hung,gm_rand] = matchMMTC(K(k),S(k),gu,dm,p_max_m,p_max_u,L,a,Bs,No);
        
        if (length(gm_hung) ~= S(k)) || (length(gm_rand) ~= S(k))
            break;
        end
        
        Em_min_hung = zeros(1,S(k));
        Em_min_rand = zeros(1,S(k));
        for s = 1:S(k)
            
            aa = gu(s);
            b = L/(Eu*Bs);
            c = -(b* log(2)* 2^(b/aa))/(aa);
            root = (-b* log(2) - aa* lambertw(-1,c))/(aa* b* log(2));
            z1 = (1/gu(s))*(2^(L/(Bs*Tu)) -1);
            
            
            if z1 < real(root) && z1 < p_max_u
                pu_opt = z1;
                pm_opt_hung = (1/gm_hung(s))*(2^(L/(Bs*Tm)) -1)*(pu_opt*gu(s) + 1);
                pm_opt_rand = (1/gm_rand(s))*(2^(L/(Bs*Tm)) -1)*(pu_opt*gu(s) + 1);
                
                Em_min_hung(s) = pm_opt_hung*Tm;
                Em_min_rand(s) = pm_opt_rand*Tm;
            else
                Em_min_hung(s) = NaN;
                Em_min_rand(s) = NaN;
                count = count +1;
                
            end
        end
        
        Em_sum_hung(i,k) = sum(Em_min_hung,'includenan');
        Em_sum_rand(i,k) = sum(Em_min_rand,'includenan');
        
    end
end
toc;

% Em_min_hung(Em_min_hung == 0) = NaN;
% Em_min_rand(Em_min_rand == 0) = NaN;
% Em_min_hung(any(isnan(Em_min_hung),2),:) = [];
% Em_min_rand(any(isnan(Em_min_rand),2),:) = [];

Em_sum_hung(any(isnan(Em_sum_hung),2),:) = [];
Em_sum_rand(any(isnan(Em_sum_rand),2),:) = [];

Em_hung_mean = mean(Em_sum_hung);
Em_rand_mean = mean(Em_sum_rand);

plot(S,sort(Em_hung_mean),'-*');
hold on
plot(S,sort(Em_rand_mean),'-*');
set(gca, 'YScale', 'log');
xlabel('Resource blocks');
ylabel('Average mMTC sum energy consumption (J)')
legend('Proposed Matching','Random Matching','Location','northwest')
hold off


