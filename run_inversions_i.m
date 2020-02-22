
% RUN_INVERSIONS_I  Optimize exponential dist. reg. wrt lambda.
% Author:           Timothy Sipkens, 2019-05-28
%=========================================================================%

if iscell(phantom.Sigma)
    Gd = phantom.Sigma{1};
else
    Gd = phantom.Sigma;
end

% [~,Gd] = phantom.p2cov(phantom.p,phantom.modes);
% Gd = Gd{2};

[x_ed_lam,lambda_ed_lam,out_ed_lam] = ...
    optimize.exp_dist_op(...
    Lb*A,Lb*b,[0.1,10],Gd,...
    grid_x.elements(:,2),grid_x.elements(:,1),x0);
disp('Process complete.');
disp(' ');

chi.ed = norm(x_ed_lam-x0);


