% function pipe_flow
clear;  clc;
dim_mat = [1 0 0 1 0;-3 1 1 -1 1;0 -1 0 -1 0];
m=5;

rho_0 = 1.2;    tol1 = 0.01;
u_0 = 10;      tol2 = 1e-4;
% u_0=0.01;
d_0 = 1.5;    tol3 = 1e-4;
mu_0 = 1e-6;    tol4 = 1e-7;
eps_0 = 2.5e-6; tol5 = 1e-6;
% eps_0 = 1e-3;   tol5 = 1e-6;

Q_nom = [rho_0,u_0,d_0,mu_0,eps_0];

tolerance = [-tol1, -tol2, -tol3, -tol4, -tol5;
    tol1, tol2, tol3, tol4, tol5];

f_nom = friction(Q_nom);
for i=1:2
    for j=1:5
        Q_input = Q_nom;
        Q_input(j) = Q_nom(j) + tolerance(i,j);
        Q_input
        f(i,j) = friction(Q_input);
    end
end
%%


f = friction(Q_nom);
del_p = 1/2 * f * Q_nom(2)^2 / (Q_nom(3) * 9.8);

disp(['f = ',num2str(f),', del_p/L = ', num2str(del_p)]);

% end
