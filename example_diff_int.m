
% clear;
% close all
% clc;


a0 = 1;

r = 1:10;

x0 = a0.*r;
b = 1/2.*a0.*(r.^2);
% epsilon = 2.*a0.*randn(1,length(r));
b_meas = b+epsilon;
x = x0+gradient(epsilon);
x1 = gradient(b_meas);

figure(1);
plot(b)
hold on;
plot(b_meas,'o-');
hold off;

figure(2);
plot(x0);
hold on;
plot(x,'o-');
plot(x1,'o-');
hold off;



