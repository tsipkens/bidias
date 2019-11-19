
% RUN_INVERSIONS_I  Multi-dimensional optimized exponential rotated and Tikhonov regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%


%-- Exponential, normal -------%
guess = [1.1,1.2,0.4]; % [lambda, sm, sd]
disp('Performing exponential distance regularization...');
[x_exp_opt,lambda_exp_opt,out_exp_opt] = optimize.exp_dist_opx(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');

chi.exp_opt = norm(x0-x_exp_opt);


%-- Exponential, normal (diff. start) -------%
guess = [1.4,0.9,0.4]; % [lambda, sm, sd]
disp('Performing exponential distance regularization...');
[x_exp_opt2,lambda_exp_opt2,out_exp_opt2] = optimize.exp_dist_opx(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');

chi.exp_opt2 = norm(x0-x_exp_opt2);


%-- Exponential rotated ------%
guess = [1.4,0.8,0.3,2.6]; % [lambda, sm, sd, Dm]
disp('Performing exponential distance regularization...');
[x_exp_rot_opt,lambda_exp_rot_opt,out_exp_rot_opt] = ...
    optimize.exp_dist_opxr(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');

chi.exp_rot_opt = norm(x0-x_exp_rot_opt);





