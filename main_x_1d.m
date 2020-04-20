
% MAIN_X_1D  Considers experiments of 1D inversion of aerosol distributions.
% Author: Timothy Sipkens, 2020-04-11
%=========================================================================%

clear; close all; clc;

disp('Reading in data...');
% [data,d_star,prop_dma] = io.import_dma('..\data\2020 UBC - Mask testing\NaCl8Apr2020.TXT');
[data,d_star,prop_dma] = io.import_dma('..\data\2020 UBC - Mask testing\2020-04-16\15Apr2020b.TXT');
disp('Complete.');
disp(' ');


%%
idx_1 = 19; % downstream
idx_2 = 20; % upstream / without filter

b1 = data(:,idx_1); % downstream
b1_max = max(b1);
b1 = b1;

b2 = data(:,idx_2); % upstream
b2_max = max(b2);


nd = 110;
d = logspace(log10(min(d_star)),log10(max(d_star)),nd)';

A = kernel.gen_1d_dma(d_star,d,prop_dma);

[~,Lb1] = tools.get_noise(b1,1,10^0.5/b1_max);
[~,Lb2] = tools.get_noise(b2,1,10^0.5/b2_max);



figure(1);
loglog(d_star,b1);
hold on;
loglog(d_star,b2,'--');
hold off;



%-- Least-squares ---------%
% x_lsq = invert.lsq(A,b);


%-- 1st order Tikhonov ----%
%-{
disp('Performing Tikhonov...');
[~,lambda_tk10,~,~] = invert1d.tikhonov_op(...
    Lb1*A,Lb1*b1,[1e-2,1e3],2,120);
disp(' ');
[~,lambda_tk20,out_tk2,~] = invert1d.tikhonov_op(...
    Lb2*A,Lb2*b2,[1e-2,1e3],2,120);

%%
lambda_tk1 = 3*lambda_tk10; % account for under-regularization for Bayes factor approach
lambda_tk2 = 3*lambda_tk20;

[x_tk1,~,~,Gpo_inv_tk1] = invert1d.tikhonov(...
    Lb1*A,Lb1*b1,lambda_tk1,2);
Gpo_tk1 = inv(Gpo_inv_tk1);

[x_tk2,~,~,Gpo_inv_tk2] = invert1d.tikhonov(...
    Lb2*A,Lb2*b2,lambda_tk2,2);
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
loglog(d,x1,'b');
hold on;
loglog(d,x1+sqrt(diag(Gpo1)),'b--');
loglog(d,max(x_tk1-sqrt(diag(Gpo1)),0),'b--');

loglog(d,x2,'r');
loglog(d,x2+sqrt(diag(Gpo2)),'r--');
loglog(d,max(x_tk2-sqrt(diag(Gpo2)),0),'r--');
hold off;
ylim([0.1*min(min(b1(b1>0)),min(b2(b2>0))),inf]);

% semilogx(d,x1,'r');
% hold on;
% semilogx(d,x1+sqrt(diag(Gpo1)),'r--');
% semilogx(d,max(x1-sqrt(diag(Gpo1)),0),'r--');
% hold off;


%-{

a1 = mvnrnd(x1,diag(diag(Gpo1)),1e4);
a2 = mvnrnd(x2,diag(diag(Gpo2)),1e4);
r = mean(a1./a2);



figure(3);
semilogx(d,100.*r,'k');
hold on;
hh = area(d,100.*(r+2.*std(a1./a2)));
hh.EdgeColor = [1,0.7,0.7];
hh.FaceColor = [1,0.9,0.9];

hl = area(d,100.*(r-2.*std(a1./a2)));
hl.EdgeColor = [1,0.7,0.7];
hl.FaceColor = [1,1,1];
semilogx(d,100.*r,'k');

semilogx(d_star,b1./b2.*100,'k.');
hold off;
ylim([0,5]);
xlim([32,max(d_star)]);

ylabel('Penetration / [%]');
xlabel('Mobility diameter / [nm]');

%}

