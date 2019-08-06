function [y] = ishigami(x1,x2,x3,a,b)
% function y = ishigami(x1,x2,x3)
% a = 7;    b = 0.1;
% a = 7;    b = 0.05;

y = sin(x1) + a* sin(x2)^2 + b*x3^4 * sin(x1);
% D = (a^2)/8 + (b*pi^4) / 5 + (b^2 * pi^8)/ 18 + 1/2;
% D1 = (b*pi^4) / 5 + (b^2 * pi^8)/50 + 1/2;
% D2 = (a^2)/8;
% D3 = 0;

end