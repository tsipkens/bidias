
% x = grid_x.reshape(x_plot);

[~,t1,t2] = grid_x.vectorize();

corr2cov = @(sigma,R) diag(sigma)*R*diag(sigma);

fun3 = @(x) x(1).*mvnpdf(log10([t1,t2]),[x(2),x(3)],...
    corr2cov(x(4).*[x(5),1],[1,x(6);x(6),1]));
x0 = [max(x_plot),0,2.3,0.3,3,0.99];
    % [C,mg,dg,sigma,Dm,corr]

x4 = lsqnonlin(@(x) fun3(x)-x_plot, x0, ...
    [0,-10,-10,0,0,-1],[inf,10,10,10,3,1]);

mu = [x4(2),x4(3)];
sigma = x4(4).*[x4(5),1];
Sigma = corr2cov(sigma,[1,x4(6);x4(6),1]);

figure(94);
t3 = fun3(x4);
grid_x.plot2d_marg(t3);

tools.plot_ellipse(mu,Sigma,1,'r')
tools.plot_ellipse(mu,Sigma,2,'r')
tools.plot_ellipse(mu,Sigma,3,'r')
% tools.plot_ellipse([x4(2),x4(3)],...
%     2.*x4(4).*[x4(5),1],x4(5),'r')
