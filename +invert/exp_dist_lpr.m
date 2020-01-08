
% EXP_DIST_LPR A helper function to 'exp_dist' to compute prior covariance.
% Author:   Timothy Sipkens, 2019-12-11
%=========================================================================%

function [Lpr,D,Gpr] = exp_dist_lpr(grid_d,m,Gd)


%-- Compute distances between elements -----------------------------------%
[vec_d1,vec_d2] = ndgrid(grid_d,grid_d);
[vec_m1,vec_m2] = ndgrid(m,m);

Gd_inv = inv(Gd);
drm = log10(vec_m1)-log10(vec_m2);
drd = log10(vec_d1)-log10(vec_d2);
D = sqrt(drm.^2.*Gd_inv(1,1)+...
    2.*drd.*drm.*Gd_inv(1,2)+...
    drd.^2.*Gd_inv(2,2)); % distance


%-- Compute prior covariance matrix --------------------------------------%
Gpr = exp(-D);

Gpr_inv = pinv(Gpr);
[Lpr,~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory

Lpr(D>1.75) = 0;
Lpr = sparse(Lpr);


end

