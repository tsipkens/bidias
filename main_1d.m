
% MAIN_1D  Considers simulations of 1D inversion of aerosol distributions.
% Author: Timothy Sipkens, 2020-04-11
%=========================================================================%

d = logspace(log10(10),log10(1e3),400)';
d_star = logspace(log10(10),log10(1e3),150)';

prop_dma = kernel.prop_dma;


A = kernel.gen_1d_dma(d_star,d,prop_dma);


mu_d = 200;
s_d = 1.25;
x0 = normpdf(log10(d),log10(mu_d),log10(s_d));
x0_alt = d.*log(10).*lognpdf(d,log(mu_d),log(s_d));
b0 = A*x0;

[b,Lb] = tools.get_noise(b0,1e5,1e-6);


figure(1);
semilogx(d_star,b,'.');
hold on;
semilogx(d_star,b0,'Color',[0.5,0.5,0.5]);
hold off;


%-- Least-squares ---------%
% x_lsq = invert.lsq(A,b);


%-- 1st order Tikhonov ----%
disp('Performing Tikhonov...');
lambda_tk1 = 2e1;
[x_tk1,~,~,Gpo_inv_tk1] = invert1d.tikhonov(Lb*A,Lb*b,lambda_tk1,1);
Gpo_tk1 = inv(Gpo_inv_tk1);
disp('Complete.');
disp(' ');


%-- Exponential distance --%
disp('Performing exponential distance...');
lambda_ed = 5e0;
ld = 1.1.*log10(s_d);
[x_ed,~,~,Gpo_inv_ed] = invert1d.exp_dist(Lb*A,Lb*b,lambda_ed,ld,d);
Gpo_ed = inv(Gpo_inv_ed);
disp('Complete.');
disp(' ');



figure(2);
semilogx(d,x0);
hold on;
semilogx(d,x_tk1);
semilogx(d,x_tk1+sqrt(diag(Gpo_tk1)),'r--');
semilogx(d,max(x_tk1-sqrt(diag(Gpo_tk1)),0),'r--');

semilogx(d,x_ed);
semilogx(d,x_ed+sqrt(diag(Gpo_ed)),'y--');
semilogx(d,max(x_ed-sqrt(diag(Gpo_ed)),0),'y--');
hold off;


