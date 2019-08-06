% hw prob 2
% sobol indices
clear;  clc;

M = 10000;      % size of sample
p = 3;          % the number of parameters

for iteration=1:10000
    
    % generate random A, B matrices using probability distribution
    A = -pi + 2*pi *rand(M,p);
    B = -pi + 2*pi *rand(M,p);
    
    C = zeros(M,p,p);
    for i=1:p
        C(:,:,i) = B;    C(:,i,i) = A(:,i);
    end
    
    y_A = zeros(M,1);   y_B = zeros(M,1); y_C = zeros(M,p);
    
    
    % for ishigami function
    a = 7;  b = 0.01;
    
    for i=1:M
        y_A(i)  = ishigami(  A(i,1),  A(i,2),  A(i,3),a,b);
        y_B(i)  = ishigami(  B(i,1),  B(i,2),  B(i,3),a,b);
        for j=1:p
            y_C(i,j) = ishigami( C(i,1,j), C(i,2,j), C(i,3,j), a, b);
        end
    end
    
    f0_sq = 1/M*sum(y_A) * 1/M*sum(y_B);
    s = zeros(p,1);     s_tot = zeros(p,1);
    %    disp('computational solution');
    for i=1:p
        s(i)= (y_C(:,i)'*y_A/M - f0_sq) / (y_A'*y_A/M - f0_sq);
        s_tot(i) = 1 - ( y_B'*y_C(:,i)/M - f0_sq ) / ( y_A'*y_A/M - f0_sq);
        %     disp(['  sobol index for x',num2str(i),' = ',num2str(s(i))]);
        %     disp(['  total effect of x',num2str(i),' = ',num2str(s_tot(i))]);
    end
    
    s_save(:,iteration) = [s; s_tot];
    
end

figure();
for i=1:3
    subplot(2,3,i);
    hist(s_save(i,:),30);
    xlabel(['s_',num2str(i)]);  ylabel('count');
    title(['histogram of s_',num2str(i)]);
    
    subplot(2,3,i+3);
    hist(s_save(i+3,:),30);
    xlabel(['sT_',num2str(i)]); ylabel('count');
    title(['histogram of sT_',num2str(i)]);
end


% 95% confidence interval :
CI = [(mean(s_save') - 1.96/sqrt(iteration)*std(s_save'))' , (mean(s_save') + 1.96/sqrt(iteration)*std(s_save'))'];

% when analytic solution is available, display it
D = (a^2)/8 + (b*pi^4) / 5 + (b^2 * pi^8)/ 18 + 1/2;
D1 = (b*pi^4) / 5 + (b^2 * pi^8)/50 + 1/2;
D2 = (a^2)/8;
D3 = 0;     D12 = 0;    D23=0;  D123 = 0;
D13 = b^2*pi^8/18 - b^2*pi^8 / 50;

disp(' ');
disp('Analytic solution');
disp(['  sobol index for x1 = ',num2str(D1/D)]);
disp(['  sobol index for x2 = ',num2str(D2/D)]);
disp(['  sobol index for x3 = ',num2str(D3/D)]);

disp(['  total effect of x1 = ',num2str((D1+D12+D13+D123)/D)]);
disp(['  total effect of x2 = ',num2str((D2+D12+D23+D123)/D)]);
disp(['  total effect of x3 = ',num2str((D3+D23+D13+D123)/D)]);


% end