
% RUN_INVERSIONS_I  Multi-dimensional optimized exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%


%-- Exponential, rotated -------%
guess = [1,1/4,log10(1.4),2.6]; % [lambda, sm, sd]
disp('Optimizing exponential distance regularization...');
[x_exp_opt,lambda_exp_opt,out_exp_opt] = optimize.exp_dist_opx(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');


disp('Parametric study of exponential distance regularization...');
[x_exp_par,lambda_exp_par,out_exp_par] = optimize.exp_dist_opbf(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');

chi.exp_opt = norm(x0-x_exp_opt);





