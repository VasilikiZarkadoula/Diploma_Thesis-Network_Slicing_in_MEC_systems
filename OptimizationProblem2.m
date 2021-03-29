clear;
clc;
rng default

%%% Constants %%%

K = 10;          % total users
M = 3;          % low latency users
N = K-M;        % low energy consumption users

R_k = 10^5 + 4*10^5*rand(K,1);  % bits
C_k = 500 + 1500*rand(K,1);     % CPU cycles 
h = 10^(-3);                    % channel gain
N_0 = 10^(-9);                  % noise variance

p_max = 2;      % max transmission power of users
B = 10^6;       % system total bandwidth
F = 6*10^9;     % total computation capacity of cloud   

E = 0.3;        % energy threshold for low latency users
E2 = 0.1;       % energy threshold for low energy users
T2 = 0.1;       % latency threshold for low energy users


%%% Objective Function %%%

fun = @(x) x(1);   

%%% Initial Variable Guess %%% 

B_m = 10^3 + 10^5*rand(M,1);
B_n = 10^3 + 10^5*rand(N,1);
F_m = 10^7 + 10^8*rand(M,1);
F_n = 10^7 + 10^8*rand(N,1);
p_m = 0.4 + 0.1*rand(M,1);
p_n = 0.1 + 0.2*rand(N,1);

x0 = [0.0001; B_m; F_m; p_m; B_n; F_n;   p_n];

%%% Variable Bounds %%%

lb = zeros(size(x0));  % all variables are positive
ub = [];

%%% Linear Equality - Inequality Constraints %%%

% Constraint 10
mat1 = eye(M);
mat2 = eye(N);
A1 = zeros(K,length(x0));
for i = 1:M
    A1(i,:) = [0 zeros(1,M) zeros(1,M) mat1(i,:) zeros(1,N) zeros(1,N) zeros(1,N)];
end
for i = 1:N
    A1(i+M,:) = [0 zeros(1,M) zeros(1,M) zeros(1,M) zeros(1,N) zeros(1,N) mat2(i,:)];
end
b1 = repmat(p_max,K,1);

% Constraint 11
A2 = [0 ones(1,M) zeros(1,M) zeros(1,M) ones(1,N) zeros(1,N) zeros(1,N)];
b2 = B;

% Constraint 12
A3 = [0 zeros(1,M) ones(1,M) zeros(1,M) zeros(1,N) ones(1,N) zeros(1,N)];
b3 = F;


A = [A1; A2; A1];
b = [b1; b2; b1];

Aeq = [];
beq = [];

% %%% Non-linear Inequality Constraints %%%
% 
% temp = constraints6to9(R_k,h,N_0,C,E,E2,T2,M,N);
% nonlcon = temp;
% 
% %%% Find optimal solution %%%
% 
% options = optimoptions(@fmincon,'MaxFunctionEvaluations',50000);
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

