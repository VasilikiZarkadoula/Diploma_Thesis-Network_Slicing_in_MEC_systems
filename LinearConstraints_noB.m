function [A,b] = LinearConstraints_noB(x0_2,M,N,K,p_max,F)

% Constraint 10
mat1 = eye(M);
mat2 = eye(N);
A1 = zeros(K,length(x0_2));
for i = 1:M
    A1(i,:) = [0 zeros(1,M) mat1(i,:) zeros(1,N) zeros(1,N)];
end
for i = 1:N
    A1(i+M,:) = [0 zeros(1,M) zeros(1,M) zeros(1,N) mat2(i,:)];
end
b1 = p_max*ones(K,1);

% % Constraint 11
% A2 = [0 ones(1,M) zeros(1,M) zeros(1,M) ones(1,N) zeros(1,N) zeros(1,N)];
% b2 = B;

% Constraint 12
A3 = [0 ones(1,M) zeros(1,M) ones(1,N) zeros(1,N)];
b3 = F;


A = [A1; A3];
b = [b1; b3];