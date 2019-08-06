% hw prob 2
% sobol indices
clear;  clc;

M = 100000;      % size of sample
p = 5;          % the number of parameters

Q_nom = [1.2 10.0 1.5 1e-6 2.5e-6];
tol = [0.01 1e-4 1e-4 1e-7 1e-6];
% tol = [0.1 0.1 0.1 0.1 0.1];


% generate random A, B matrices using probability distribution
% assum the parameter is uniformly distributed
for i=1:M 
    A(i,:) = Q_nom-tol + 2*tol.*rand(1,p);
    B(i,:) = Q_nom-tol + 2*tol.*rand(1,p);
end

% C = zeros(M,p,p);
C1 = B;     C1(:,1) = A(:,1);
C2 = B;     C2(:,2) = A(:,2);
C3 = B;     C3(:,3) = A(:,3);
C4 = B;     C4(:,4) = A(:,4);
C5 = B;     C5(:,5) = A(:,5);

y_A = zeros(M,1);   
y_B = zeros(M,1); 
y_C1 = zeros(M,1);
y_C2 = zeros(M,1);
y_C3 = zeros(M,1);
y_C4 = zeros(M,1);
y_C5 = zeros(M,1);


for i=1:M
    y_A(i)  = 0.5*A(i,2)^2 / (A(i,3)*9.8)*friction(A(i,:));
    y_B(i)  = 0.5*B(i,2)^2 / (B(i,3)*9.8)*friction(B(i,:));
    
    y_C1(i) = 0.5*C1(i,2)^2 / (C1(i,3)*9.8)*friction(C1(i,:));
    y_C2(i) = 0.5*C2(i,2)^2 / (C2(i,3)*9.8)*friction(C2(i,:));
    y_C3(i) = 0.5*C3(i,2)^2 / (C3(i,3)*9.8)*friction(C3(i,:));
    y_C4(i) = 0.5*C4(i,2)^2 / (C4(i,3)*9.8)*friction(C4(i,:));
    y_C5(i) = 0.5*C5(i,2)^2 / (C5(i,3)*9.8)*friction(C5(i,:));
    
end

f0_sq = 1/M*sum(y_A) * 1/M*sum(y_B);
s = zeros(p,1);     s_tot = zeros(p,1);


s(1) = (y_C1'*y_A/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s(2) = (y_C2'*y_A/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s(3) = (y_C3'*y_A/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s(4) = (y_C4'*y_A/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s(5) = (y_C5'*y_A/M - f0_sq) / (y_A'*y_A/M - f0_sq);

s_tot(1) = 1 - (y_C1'*y_B/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s_tot(2) = 1 - (y_C2'*y_B/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s_tot(3) = 1 - (y_C3'*y_B/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s_tot(4) = 1 - (y_C4'*y_B/M - f0_sq) / (y_A'*y_A/M - f0_sq);
s_tot(5) = 1 - (y_C5'*y_B/M - f0_sq) / (y_A'*y_A/M - f0_sq);



s'
s_tot'





% end