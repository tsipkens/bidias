
% RUN_INVERSIONS_I  Optimize exponential distance regularization w.r.t. lambda.
% Author: Timothy Sipkens, 2019-05-28
%=========================================================================%

if iscell(phantom.Sigma)
    Gd = phantom.Sigma{1};
elseif isempty(phantom.Sigma)
    p = phantom.p; ll = 2;
    Gd = inv([(1/log10(p(ll).smd))^2,...
        -p(ll).Dm/log10(p(ll).smd)^2;...
        -p(ll).Dm/log10(p(ll).smd)^2,...
        1/log10(p(ll).sg)^2+p(ll).Dm^2/log10(p(ll).smd)^2]);
else
    Gd = phantom.Sigma;
end

%-- Gd properties -----------------%
l1 = sqrt(Gd(1,1));
l2 = sqrt(Gd(2,2));
R12 = Gd(1,2)/(l1*l2);
Dm = Gd(1,2)/Gd(2,2); % s1*R12/s2
%----------------------------------%

[x_ed_lam,lambda_ed_lam,out_ed_lam] = ...
    optimize.exp_dist_op(...
    Lb*A,Lb*b,[0.1,10],Gd,...
    grid_x,[],x0);
disp('Process complete.');
disp(' ');

eps.ed = norm(x_ed_lam-x0);


