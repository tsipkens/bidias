
% MAIN_X_1D  Considers experiments of 1D inversion of aerosol distributions.
% Author: Timothy Sipkens, 2020-04-11
%=========================================================================%

[data,d_star,prop_dma] = ...
    io.import_1dma('..\data\2020 UBC - Mask testing\NaCl8Apr2020.TXT');
prop_dma.G_DMA = 3.0256;
prop_dma.bet = 0.1000;
prop_dma.del = 0;
prop_dma.L = prop_dma.L./100;
prop_dma.p = prop_dma.p./101.325;
prop_dma.R1 = prop_dma.R1/100;
prop_dma.R2 = prop_dma.R2/100;


%%
b1 = data(:,14); % downstream
b1_max = max(b1);
b1 = b1;

b2 = data(:,11); % upstream
b2_max = max(b2);
b2 = b2;


nd = 400;
d = logspace(log10(10),log10(1e3),nd)';

A = kernel.gen_1d_dma(d_star,d,prop_dma);


[~,Lb1] = tools.get_noise(b1,1,1e-1);
[~,Lb2] = tools.get_noise(b2,1,1e-1);



figure(1);
semilogx(d_star,b1);
hold on;
semilogx(d_star,b2,'--');
hold off;



%-- Least-squares ---------%
% x_lsq = invert.lsq(A,b);


%-- 1st order Tikhonov ----%
%-{
disp('Performing Tikhonov...');
lambda_tk1 = 1.5e0;
[x_tk1,~,~,Gpo_inv_tk1] = invert1d.tikhonov(Lb1*A,Lb1*b1,lambda_tk1,1);
Gpo_tk1 = inv(Gpo_inv_tk1);

lambda_tk2 = lambda_tk1;
[x_tk2,~,~,Gpo_inv_tk2] = invert1d.tikhonov(Lb2*A,Lb2*b2,lambda_tk2,1);
Gpo_tk2 = inv(Gpo_inv_tk2);
disp('Complete.');
disp(' ');

x1 = x_tk1;
x2 = x_tk2;
Gpo1 = Gpo_tk1;
Gpo2 = Gpo_tk2;
%}


%-- Exponential distance --%
%{
disp('Performing exponential distance...');
lambda_ed1 = 1e0;
ld1 = 0.05;
[x_ed1,~,~,Gpo_inv_ed1] = invert1d.exp_dist(Lb1*A,Lb1*b1,lambda_ed1,ld1,d);
Gpo_ed1 = inv(Gpo_inv_ed1);

lambda_ed2 = 1e0;
ld2 = 0.05;
[x_ed2,~,~,Gpo_inv_ed2] = invert1d.exp_dist(Lb2*A,Lb2*b2,lambda_ed2,ld2,d);
Gpo_ed2 = inv(Gpo_inv_ed2);
disp('Complete.');
disp(' ');

x1 = x_ed1;
x2 = x_ed2;
Gpo1 = Gpo_ed1;
Gpo2 = Gpo_ed2;
%}



figure(2);
% semilogx(d,x_tk1,'b');
% hold on;
% semilogx(d,x_tk1+sqrt(diag(Gpo_tk1)),'b--');
% semilogx(d,max(x_tk1-sqrt(diag(Gpo_tk1)),0),'b--');

semilogx(d,x1,'r');
hold on;
semilogx(d,x1+sqrt(diag(Gpo1)),'r--');
semilogx(d,max(x1-sqrt(diag(Gpo1)),0),'r--');
hold off;


%-{

a1 = mvnrnd(x1,Gpo_ed1,1e4);
a2 = mvnrnd(x2,Gpo_ed2,1e4);
r = mean(a1./a2);


figure(3);
semilogx(d,r,'k');
hold on;
hh = area(d,r+2.*std(a1./a2));
hh.EdgeColor = [1,0.7,0.7];
hh.FaceColor = [1,0.9,0.9];

hl = area(d,r-2.*std(a1./a2));
hl.EdgeColor = [1,0.7,0.7];
hl.FaceColor = [1,1,1];
semilogx(d,r,'k');
hold off;
ylim([0,0.3]);
xlim([30,800]);

%}

