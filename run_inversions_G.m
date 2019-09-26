
% RUN_INVERSIONS_G  Optimize exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_Tk1,lambda_Tk1,out_Tk1] = invert.optimize_tikhonov(...
    Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,1);
t.Tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk1 = norm(x0-x_Tk1);


%-- Exponential rotated ------%
s1 = 2.0;
s2 = 0.4;
dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
theta = -atan2(1,2.5);
Lex = diag([1/s1,1/s2])*...
    [cos(theta),-sin(theta);sin(theta),cos(theta)];
% lambda_expRot = 1.8;

disp('Performing rotated exponential distance regularization...');
[x_expRot,lambda_expRot,out_expRot] = invert.optimize_exp_dist(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    [1e-2,1e2],x0,Lex);
disp('Inversion complete.');
disp(' ');

chi.exp_dist = norm(x0-x_expRot);


%-- Calculate Bayes factor ----------%
Fi_pr = []; % fit
Fi_b = [];
Gpo_inv_logdet = [];
Gpo_logdet = [];
n = prod(n_x);
tools.textbar(0);
for ii=1:length(out_Tk1)
    Gpo_logdet(ii) = tools.logdet(inv(out_Tk1(ii).Gpo_inv));
    Fi_pr(ii) = -1/2.*(norm(out_Tk1(ii).Lpr*(out_Tk1(ii).x))^2);
    Fi_b(ii) = -1/2.*(norm(Lb*(A*out_Tk1(ii).x-b))^2);
    tools.textbar(ii/length(out_Tk1));
end
Ci = 1/2.*(Gpo_logdet+...
    2*n.*log([out_Tk1.lambda]));
Fi = Fi_pr+Fi_b;
Bi = Fi+Ci;
% semilogx([out_Tk1.lambda],[Ci;Fi;Bi]');

