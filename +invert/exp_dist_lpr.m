
% EXP_DIST_LPR A helper function to 'exp_dist' to compute prior covariance.
% Author:   Timothy Sipkens, 2019-12-11
%=========================================================================%

function [Lpr] = exp_dist_lpr(d_vec,m_vec,lambda,Gd)


%-- Compute distances between elements -----------------------------------%
[vec_d1,vec_d2] = ndgrid(d_vec,d_vec);
[vec_m1,vec_m2] = ndgrid(m_vec,m_vec);

Gd_inv = inv(Gd);
drm = log10(vec_m1)-log10(vec_m2);
drd = log10(vec_d1)-log10(vec_d2);
d = sqrt(drm.^2.*Gd_inv(1,1)+...
    2.*drd.*drm.*Gd_inv(1,2)+...
    drd.^2.*Gd_inv(2,2)); % distance


%-- Compute prior covariance matrix --------------------------------------%
Gpr = exp(-d);

Gpr_inv = pinv(Gpr);
[Lpr,~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory
Lpr = lambda.*Lpr./max(max(Lpr));

Lpr(d>1.5) = 0;
Lpr = sparse(Lpr);


end

