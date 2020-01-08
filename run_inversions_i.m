
% RUN_INVERSIONS_I  Multi-dimensional optimized exponential rotated and exp. dist. regularization
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

% if iscell(phantom.Sigma)
%     Gd = phantom.Sigma{1};
% else
%     Gd = phantom.Sigma;
% end

[~,Gd] = phantom.p2cov(phantom.p,phantom.modes);
Gd = Gd{2};

[x_ed_lam,lambda_ed_lam,out_ed_lam] = ...
    optimize.exp_dist_op(Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    [0.1,10],x0,Gd);
disp('Process complete.');
disp(' ');

chi.ed = norm(x_ed_lam-x0);


%{
guess = [1.3,1/4,log10(1.8),0.84]; % [lambda, ratio, ld, corr]
disp('Optimizing exponential distance regularization (least-sqaures)...');
[x_ed_opt,lambda_ed_opt,out_ed_opt] = optimize.exp_dist_opx(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    guess,x0); 
disp('Inversion complete.');
disp(' ');
%}

%{
disp('Parametric study of exponential distance regularization (brute force)...');
[x_ed_par,lambda_ed_par,out_ed_par] = optimize.exp_dist_opbf(...
    Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    x0); 
disp('Inversion complete.');
disp(' ');
%}


%{
if iscell(phantom.Sigma)
    Gd = phantom.Sigma{1};
else
    Gd = phantom.Sigma;
end
[out_ed_corr] = ...
    optimize.exp_dist_op1d(Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
    lambda_ed_lam,x0,Gd);
%}

%{
%-- Zeroth-order Tikhonov regularization --%
%   Lower limit for correlation lengths.
[x_tk0,D_tk0,L_tk0,Gpo_tk0] = invert.tikhonov(Lb*A,Lb*b,n_x(1),...
    lambda_ed_lam,0);
%}
