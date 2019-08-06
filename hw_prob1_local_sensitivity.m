%% UQ sensitivity homework
% problem 1
% local sensitivity - adjoint method

x = [-0.5 0.5 1]';
B = [-0.1 0.5 0.0;-0.5 -0.1 0.3;0.0 -0.3 -0.2];
A_inv = expm(B);

t = 0.0:0.1:20;
y = zeros(3,length(t));
for i=1:length(t)
    y(:,i) = expm(B*t(i))*x;
end

figure();
plot(t,y(1,:)); hold on
plot(t,y(2,:)); 
plot(t,y(3,:));

f = sqrt(sum(y(:,end).^2));


%%
x_range = -1:0.1:1;
t_f = 20;
f = zeros(length(x_range),length(x_range),length(x_range));

for i=1:length(x_range)
    for j=1:length(x_range)
        for k=1:length(x_range)
            x=[x_range(i), x_range(j), x_range(k)];
            y_f = expm(B*t_f)*x';
            f(i,j,k) = sqrt(sum(y_f.^2));
        end
    end
end
        

%% UQ sensitivity homework
% problem 1
% local sensitivity - adjoint method
clear; clc;

B = [-0.1 0.5 0.0;-0.5 -0.1 0.3;0.0 -0.3 -0.2];
x_range = -1:0.01:1;    % the range of one x component
t_f = 20;               % final time 
y_f = zeros(3,length(x_range));
df_dy = zeros(3,length(x_range));
f   = zeros(1,length(x_range));

for i=1:length(x_range)
    x(1) = -0.1;  % discretize only one element of x
    x(2) = 0.5;         % while two others are fixed
    x(3) = x_range(i);
    
    y_f(:,i) = expm(B*t_f)*x';  % calculate y 
    f(i) = sqrt(sum(y_f(:,i).^2));
    
    % compute df/dy from analytic derivation
    df_dy(:,i) = y_f(:,i) * f(i)^(-1);
    df_dx(:,i) = df_dy(:,i)'*expm(B*t_f);
end

% numerical computation of df/dx
df_dx2 = zeros(1,length(x_range));
df_dx2(2:end-1) = (f(3:end) - f(1:end-2)) / (x_range(3)-x_range(1));
df_dx2(1) = (f(2) - f(1)) / (x_range(2) - x_range(1));
df_dx2(end) = (f(end) - f(end-1)) / (x_range(end) - x_range(end-1));


figure();
subplot(1,2,1); plot(x_range, df_dx(3,:));
xlabel('x(3)');    ylabel('df/dx');
title('analytic, x(1) = -0.1, x(2) = 0.5');
hold on
subplot(1,2,2); plot(x_range, df_dx2);
xlabel('x(3)');    ylabel('df/dx');
title('computational, x(1) = -0.1, x(2) = 0.5');



