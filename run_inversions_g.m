
% RUN_INVERSIONS_G  Optimize exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

%-- Tikhonov (1st order) -----%
disp('Performing Tikhonov (1st) regularization...');
tic;
[x_tk1,lambda_tk1,out_tk1] = optimize.tikhonov_op(...
    Lb*A,Lb*b,n_x(1),[1e-2,1e1],x0,1);
t.tk1 = toc;
disp('Inversion complete.');
disp(' ');

chi.tk1_nn = norm(x0-x_tk1);




