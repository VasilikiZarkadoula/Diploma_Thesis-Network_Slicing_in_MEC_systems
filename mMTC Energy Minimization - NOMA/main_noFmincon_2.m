close all;
clear;
clc;

K = 10;       
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
                  
Tm = 0.1:0.1:1.5;
Eu = 5e-2;
Tu = 0.01;
b = L/(Eu*Bs);

iterations = 1;
count = 0;

tic;

Em_min_hung = zeros(1,S);
Em_min_rand = zeros(1,S);
Em_min_exh = zeros(1,S);
Em_sum_hung = zeros(iterations,length(Tm));
Em_sum_rand = zeros(iterations,length(Tm));
Em_sum_exh = zeros(iterations,length(Tm));
for i = 1:iterations
    
    [dm,~,~] = distance_calc(dmin_k, dmin_b, radius_m , K/2);
    [du,~,~] = distance_calc(dmin_k, radius_m, radius_u , K/2);
    
    gu = matchURLLC(K,S,du,a,Bs,No);
    [gm_hung,gm_rand,gm_exh] = matchMMTC(K,S,gu,dm,p_max_m,p_max_u,L,a,Bs,No);
    
    if (length(gm_hung) ~= S) || (length(gm_rand) ~= S) || (length(gm_exh) ~= S)
        break;
    end
    for t = 1:length(Tm)
        for s = 1:S
            
            aa = gu(s);
            b = L/(Eu*Bs);
            c = -(b* log(2)* 2^(b/aa))/(aa);
            root = (-b* log(2) - aa* lambertw(-1,c))/(aa* b* log(2));
            z1 = (1/gu(s))*(2^(L/(Bs*Tu)) -1);
            
            
            if z1 < real(root) && z1 < p_max_u
                pu_opt = z1;
                pm_opt_hung = (1/gm_hung(s))*(2^(L/(Bs*Tm(t))) -1)*(pu_opt*gu(s) + 1);
                pm_opt_rand = (1/gm_rand(s))*(2^(L/(Bs*Tm(t))) -1)*(pu_opt*gu(s) + 1);
                pm_opt_exh = (1/gm_exh(s))*(2^(L/(Bs*Tm(t))) -1)*(pu_opt*gu(s) + 1);
                
                Em_min_hung(s) = pm_opt_hung*Tm(t);
                Em_min_rand(s) = pm_opt_rand*Tm(t);
                Em_min_exh(s) = pm_opt_exh*Tm(t);
            else
                Em_min_hung(s) = NaN;
                Em_min_rand(s) = NaN;
                Em_min_exh(s) = NaN;
                count = count +1;
                
            end
        end
        Em_sum_hung(i,t) = sum(Em_min_hung,'includenan');
        Em_sum_rand(i,t) = sum(Em_min_rand,'includenan');
        Em_sum_exh(i,t) = sum(Em_min_exh,'includenan');
    end
end
toc;

% Em_min_hung(Em_min_hung == 0) = NaN;
% Em_min_rand(Em_min_rand == 0) = NaN;
% Em_min_hung(any(isnan(Em_min_hung),2),:) = [];
% Em_min_rand(any(isnan(Em_min_rand),2),:) = [];

Em_sum_hung(any(isnan(Em_sum_hung),2),:) = [];
Em_sum_rand(any(isnan(Em_sum_rand),2),:) = [];
Em_sum_exh(any(isnan(Em_sum_exh),2),:) = [];


Em_hung_mean = mean(Em_sum_hung);
Em_rand_mean = mean(Em_sum_rand);
Em_exh_mean = mean(Em_sum_exh);

% plot(Tm,Em_hung_mean*1e+3,'k','LineWidth',1.7);
% hold on
% plot(Tm,Em_rand_mean*1e+3);
% hold on
% plot(Tm,Em_exh_mean*1e+3,'-*');

semilogy(Tm,Em_hung_mean*1e+3,'k','LineWidth',1.7);
hold on
semilogy(Tm,Em_rand_mean*1e+3);
hold on
semilogy(Tm,Em_exh_mean*1e+3,'-*')
% set(gca, 'YScale', 'log');
xlabel('mMTC latency threshold (sec)');
ylabel('Average mMTC sum energy consumption (mJ)')
legend('Proposed Matching','Random Matching','Exhaustive Search','Location','best')
hold off


