
% RUN_INVERSIONS_G  Optimize exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

% Tikhonov (1st order)
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_Tk1,lambda_Tk1,out_Tk1] = invert.optimize_tikhonov(...
    Lb*A,Lb*b,n_x(1),[1e-2,1e2],x0,1);
t.Tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.Tk1 = norm(x0-x_Tk1);


% Exponential rotated
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
