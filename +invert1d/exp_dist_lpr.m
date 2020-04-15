
% EXP_DIST_LPR  A helper function for 'exp_dist' to compute prior covariance.
% Author: Timothy Sipkens, 2019-12-11
%=========================================================================%

function [Lpr,D,Gpr] = exp_dist_lpr(ld,vec)


%-- Compute distances between elements -----------------------------------%
[vec_a,vec_b] = ndgrid(vec,vec);

dr = log10(vec_a)-log10(vec_b);
D = abs(dr./ld); % distance, abs replaces sqrt(d.^2)


%-- Compute prior covariance matrix --------------------------------------%
Gpr = exp(-D);

Gpr_inv = pinv(Gpr);
[Lpr,~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory

Lpr(D>1.75) = 0;
Lpr = sparse(Lpr);


end

