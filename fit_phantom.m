
% [x_Tk1,D_Tk1,L_Tk1] = tikhonov(Lb*A,Lb*b,n_x(1),2e-3,1);

rho = 700;
phi0(2) = 2.4; % Dm
phi0(3) = 10^2.1; % dg
phi0(1) = rho*pi/6*phi0(3)^(3-phi0(2)); % k
phi0(4) = 1.6; % sg
phi0(5) = 1.522; % sm
phi0(6) = 1;
phi0(6) = max(x_Tk1)/max(fun(phi0,grid_x)); % scale parameter

d = 0;
phi1 = fminsearch(@(phi) norm(fun(phi,grid_x)-x_Tk1)^2, phi0);
rho_fit = 6*phi1(1)/pi*phi1(3)^(phi1(2)-3);
mg_fit  = phi1(1)*phi1(3)^phi1(2);
x_fit = fun(phi1,grid_x);

phi1(1) = 9400;
phi1(2) = 2.3;
phi1(3) = 125;


figure(1);
colormap(gcf,cm);
grid_x.plot2d(x_fit);


function x_fit = fun(phi0,grid_x)

param = [];
param(1).k = phi0(1);
param(1).Dm = phi0(2);
param(1).dg = phi0(3);
param(1).sg = phi0(4);
param(1).sm = phi0(5);

[x_fit,grid_t] = gen_phantom(param,grid_x.span);
x_fit = grid_x.project(grid_t.edges,x_fit);
x_fit = phi0(6).*x_fit;

end

