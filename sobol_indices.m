% function sobol_indices()
% 
%   x* : input parameters
%   i* : iteration for x*
%   *  : the number of parameters
%   f  : function values ( dimension : length(x1) x length(x2) * ...
%                                   * length(x*))

x1 = -pi:0.05:pi;   x2 = -pi:0.05:pi;   x3 = -pi:0.05:pi;
a = 1;  b = 1;

f = zeros(length(x1),length(x2),length(x3));
disp(['size of f = ',num2str(size(f))])


for i1 = 1:length(x1)
    for i2 = 1:length(x2)
        for i3 = 1:length(x3)
            f(i1,i2,i3) = ishigami(x1(i1),x2(i2),x3(i3),a,b);
%             f(i1,i2,i3) = function you want
        end
    end    
end


VAR = mean(mean(mean(f.^2))) - mean(mean(mean(f))).^2;
VAR_1 = var(mean(mean(f,2),3));
VAR_2 = var(mean(mean(f,3),1));
VAR_3 = var(mean(mean(f,1),2));

% calculate confidence interval
mu_hat = mean(mean(mean(f)));
CI = [mu_hat + norminv(0.025)*sqrt(VAR), mu_hat + norminv(0.975)*sqrt(VAR)];


% computational solution of sobol indices
disp('computational solution');
disp(['  total variance = ',num2str(VAR)]);
disp(['  sobol index for x1 = ',num2str(VAR_1/VAR)]);
disp(['  sobol index for x2 = ',num2str(VAR_2/VAR)]);
disp(['  sobol index for x3 = ',num2str(VAR_3/VAR)]);


% when analytic solution is available, display it
D = (a^2)/8 + (b*pi^4) / 5 + (b^2 * pi^8)/ 18 + 1/2;
D1 = (b*pi^4) / 5 + (b^2 * pi^8)/50 + 1/2;
D2 = (a^2)/8;
D3 = 0;

disp('Analytic solution');
disp(['  total variance = ',num2str(D)]);
disp(['  sobol index for x1 = ',num2str(D1/D)]);
disp(['  sobol index for x2 = ',num2str(D2/D)]);
disp(['  sobol index for x3 = ',num2str(D3/D)]);


% end