
% RUN_INVERSIONS_I  Multi-dimensional optimized exponential rotated and exp. dist. regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

if iscell(phantom.Sigma)
    Gd = phantom.Sigma{1};
else
    Gd = phantom.Sigma;
end
[x_ed_lam,lambda_ed_lam,out_ed_lam] = ...
    optimize.exp_dist_op(Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    [0.1,10],x0,Gd);
disp('Process complete.');
disp(' ');

guess = [1.3,1/4,log10(1.8),0.84]; % [lambda, ratio, ld, corr]
disp('Optimizing exponential distance regularization (least-sqaures)...');
[x_ed_opt,lambda_ed_opt,out_ed_opt] = optimize.exp_dist_opx(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');


disp('Parametric study of exponential distance regularization (brute force)...');
[x_ed_par,lambda_ed_par,out_ed_par] = optimize.exp_dist_opbf(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');

chi.exp_dist = norm(x_ed_opt-x0);

% x_ed = invert.exp_dist(...
%     Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
%     2,Gd); 



