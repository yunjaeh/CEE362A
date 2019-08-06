% hw prob 2
% sobol indices
clear;  clc;


M = 100000; % size of sample
p = 3;


for iteration = 1:1000
    

A = -pi + 2*pi *rand(M,p);
B = -pi + 2*pi *rand(M,p);

C1 = B; C1(:,1) = A(:,1);
C2 = B; C2(:,2) = A(:,2);
C3 = B; C3(:,3) = A(:,3);


y_A  = zeros(M,1);
y_B  = zeros(M,1); 
y_C1 = zeros(M,1); 
y_C2 = zeros(M,1); 
y_C3 = zeros(M,1); 

% for ishigami function
a = 7;  b = 0.01;
% a=1;    b=1;

for i=1:M
    y_A(i)  = ishigami(  A(i,1),  A(i,2),  A(i,3),a,b);
    y_B(i)  = ishigami(  B(i,1),  B(i,2),  B(i,3),a,b);
    y_C1(i) = ishigami( C1(i,1), C1(i,2), C1(i,3),a,b);
    y_C2(i) = ishigami( C2(i,1), C2(i,2), C2(i,3),a,b);
    y_C3(i) = ishigami( C3(i,1), C3(i,2), C3(i,3),a,b);
end

f_02 = 1/M*sum(y_A) * 1/M*sum(y_B);

% calculate confidence interval
VAR = (y_A'*y_A/M - f_02);
s1 = (y_C1'*y_A/M - f_02) / (y_A'*y_A/M - f_02);
s2 = (y_C2'*y_A/M - f_02) / (y_A'*y_A/M - f_02);
s3 = (y_C3'*y_A/M - f_02) / (y_A'*y_A/M - f_02);
s_tot1 = 1 - (y_B'*y_C1/M - f_02)/(y_A'*y_A/M - f_02);
s_tot2 = 1 - (y_B'*y_C2/M - f_02)/(y_A'*y_A/M - f_02);
s_tot3 = 1 - (y_B'*y_C3/M - f_02)/(y_A'*y_A/M - f_02);



% computational solution of sobol indices
% disp('computational solution');
% disp(['  total variance = ',num2str(VAR)]);
% disp(['  sobol index for x1 = ',num2str(s1)]);
% disp(['  sobol index for x2 = ',num2str(s2)]);
% disp(['  sobol index for x3 = ',num2str(s3)]);

s(iteration,1) = s1;
s(iteration,2) = s2;
s(iteration,3) = s3;
s(iteration,4) = s_tot1;
s(iteration,5) = s_tot2;
s(iteration,6) = s_tot3;

end

figure();
subplot(2,3,1);
hist(s(:,1),30);
subplot(2,3,2);
hist(s(:,2),30);
subplot(2,3,3);
hist(s(:,3),30);

subplot(2,3,4);
hist(s(:,4),30);
subplot(2,3,5);
hist(s(:,5),30);
subplot(2,3,6);
hist(s(:,6),30);

% when analytic solution is available, display it
D = (a^2)/8 + (b*pi^4)/5 + (b^2 * pi^8)/ 18 + 1/2;
D1 = (b*pi^4)/5 + (b^2*pi^8)/50 + 1/2;
D2 = (a^2)/8;
D3 = 0;
D13 = (b^2 * pi^8)/ 18 - (b^2 * pi^8) /50;


disp('Analytic solution');
disp(['  total variance = ',num2str(D)]);
disp(['  sobol index for x1 = ',num2str(D1/D)]);
disp(['  sobol index for x2 = ',num2str(D2/D)]);
disp(['  sobol index for x3 = ',num2str(D3/D)]);
disp(' ');
disp(['  total effect of x1 = ',num2str((D1+D13)/D)]);
disp(['  total effect of x2 = ',num2str(D2/D)]);
disp(['  total effect of x3 = ',num2str((D3+D13)/D)]);
mean(s)

% end