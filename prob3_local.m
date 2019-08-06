% problem 3
% local sensitivity
clear; clc; close all;
rho_0 = 1.2;
mu_0 = 1e-6;
eps_0 = 2.5e-6;

% % design 1
% u_0 = 10.0;   d_0 = 1.5;

% % design 2
% u_0 = 1.0;    d_0 = 0.1;

% % design 3
%u_0 = 0.1;      d_0 = 0.001;

Q_nom = [rho_0 u_0 d_0 mu_0 eps_0];
tol = [0.01 1e-4 1e-4 1e-7 1e-6];

% % part (ii)
eps_0 = 1.5e-4;
Q_nom(5) = eps_0;
tol(5) = 5e-5;



f_nom = 1/2*friction(Q_nom)*u_0^2 / (d_0 * 9.8); 
N = 100;    % the length of parameter

x_save = zeros(5,N);
f_save = zeros(5,N);
df_dx_save = zeros(5,N);
df_dx = zeros(1,N);

figure();
for param=1:5
    x = linspace(Q_nom(param)-tol(param),Q_nom(param)+tol(param),N);
    for i=1:N
        Q_input = Q_nom;
        Q_input(param) = x(i);
        f(i) = 1/2*friction(Q_input)*Q_input(2)^2 / (Q_input(3) * 9.8); 
    end
    df_dx(2:end-1) = (f(3:end) - f(1:end-2))/(x(3)-x(1));
    df_dx(1)       = (f(2) - f(1))/(x(2)-x(1));
    df_dx(end)     = (f(end) - f(end-1))/(x(end)-x(end-1));
    
    x_save(param,:)     = x;
    f_save(param,:)     = f;
    df_dx_save(param,:) = df_dx;
   
    subplot(2,5,param);
    plot(x,f);
    xlim([min(x) max(x)]);  grid on
    xlabel('x');    ylabel('f (QoI)');
    
    subplot(2,5,param+5);
    plot(x,df_dx);
    xlim([min(x) max(x)]); grid on
    xlabel('x');    ylabel('df/dx');
end

max(f_save')
min(f_save')

var_max = (max(f_save') - f_nom)./f_nom*100
var_min = (min(f_save') - f_nom)./f_nom*100


subplot(2,5,1);     xlabel('\rho_F');
subplot(2,5,6);     xlabel('\rho_F');

subplot(2,5,2);     xlabel('U_F');
subplot(2,5,7);     xlabel('U_F');

subplot(2,5,3);     xlabel('D_P');
subplot(2,5,8);     xlabel('D_P');

subplot(2,5,4);     xlabel('\mu_F');
subplot(2,5,9);     xlabel('\mu_F');

subplot(2,5,5);     xlabel('\epsilon_P');
subplot(2,5,10);    xlabel('\epsilon_P');

% -0.0029   -0.0037    0.7926    0.1973   -0.0037
%   0.0997    0.0979    0.8279    0.2767    0.0978
%   
%     0.0072    0.0071    0.7913    0.2123    0.0072
%    -0.0003   -0.0032    0.8008    0.1889   -0.0033