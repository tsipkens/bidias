
% RUN_INVERSIONS_G  Optimize exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(...
    Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,1);
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1 = norm(x0-x_tk1);


%-- Exponential, normal ------%
s1 = 0.3;
s2 = 0.1;
Gd = diag([s1^2,s2^2]);

disp('Performing exponential distance regularization...');
[x_exp,lambda_exp,out_exp] = optimize.exp_dist_op(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    [1e-2,1e2],x0,Gd);close 
disp('Inversion complete.');
disp(' ');

chi.exp = norm(x0-x_exp);


%-- Exponential rotated ------%
s1 = 0.8;
s2 = 0.3;
Dm = 1.92; % used to specify correlation
R12 = 0.99;%Dm*s2/s1; % correlation
Gd2 = [s1^2,R12*(s1*s2);R12*(s1*s2),s2^2];

disp('Performing rotated exponential distance regularization...');
[x_exp_rot,lambda_exp_rot,out_exp_rot] = optimize.exp_dist_op(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    [1e-2,1e2],x0,Gd2);
disp('Inversion complete.');
disp(' ');

chi.exp_rot = norm(x0-x_exp_rot);


%-- Calculate Bayes factor ----------%
%{
Fi_pr = []; % fit
Fi_b = [];
Gpo_inv_logdet = [];
Gpo_logdet = [];
n = prod(n_x);
tools.textbar(0);
for ii=1:length(out_tk1)
    Gpo_logdet(ii) = tools.logdet(inv(out_exp_rot(ii).Gpo_inv));
    Fi_pr(ii) = -1/2.*(norm(out_exp_rot(ii).Lpr*out_exp_rot(ii).x)^2);
    Fi_b(ii) = -1/2.*(norm(Lb*(A*out_exp_rot(ii).x-b))^2);
    Axb(ii) = -1/2.*(norm(A*out_exp_rot(ii).x-b)^2);
    tools.textbar(ii/length(out_exp_rot));
end
Ci = 1/2.*(Gpo_logdet+...
    2*n.*log([out_exp_rot.lambda]));
Fi = Fi_pr+Fi_b;
Bi = Fi+Ci;
semilogx([out_exp_rot.lambda],[Ci;Fi;Bi]');
%}



