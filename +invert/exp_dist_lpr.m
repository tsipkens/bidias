
% EXP_DIST_LPR  A helper function for 'exp_dist' to compute prior covariance.
% Author: Timothy Sipkens, 2019-12-11
%=========================================================================%

function [Lpr,D,Gpr] = exp_dist_lpr(Gd,grid_vec2,vec1)

%-- Parse inputs ---------------------------------------------------------%
if isa(grid_vec2,'Grid') % if a grid is supplied (`vec1` is unused)
    vec2 = grid_vec2.elements(:,2);
    vec1 = grid_vec2.elements(:,1);
else % if vectors of elements are supplied (`vec1` is used unchanged)
    vec2 = grid_vec2;
end
%-------------------------------------------------------------------------%


%-- Compute distances between elements -----------------------------------%
[vec2_a,vec2_b] = ndgrid(vec2,vec2); % for differences in 2nd dim
[vec1_a,vec1_b] = ndgrid(vec1,vec1); % for differences in 1st dim

Gd_inv = inv(Gd);
dr1 = log10(vec1_a)-log10(vec1_b);
dr2 = log10(vec2_a)-log10(vec2_b);
D = sqrt(dr1.^2.*Gd_inv(1,1)+...
    2.*dr2.*dr1.*Gd_inv(1,2)+...
    dr2.^2.*Gd_inv(2,2)); % distance


%-- Compute prior covariance matrix --------------------------------------%
Gpr = exp(-D);

Gpr_inv = pinv(Gpr);
[Lpr,~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory

Lpr(D>1.75) = 0;
Lpr = sparse(Lpr);


end

