
% EXP_DIST  Regularization based on the exponential of the distance between elements/pixels.
% Author: Timothy Sipkens, 2018-10-22
% 
% Inputs:
%   A           Model matrix or Lb*A
%   b           Data or Lb*b
%   lambda      Regularization parameter
%   Gd          Mass-mobility covariance matrix, used to calculate
%               Mahalanobis distance (Optional, default: identity matrix)
%   grid_vec2	Position of elements in mobility space, grid.elements(:,2)
%   vec1        Position of elements in mass space, grid.elements(:,1)
%   xi          Initial guess for solver (Optional, default: empty)
%   solver      Least-squares solver to use (e.g. 'interior-point', see 'invert.lsq' function)
%
% Outputs:
%   x           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   Lpr0        Cholesky factorization of prior covariance
%   Gpo_inv     Inverse of the posterior covariance
%=========================================================================%

function [x,D,Lpr0,Gpo_inv] = exp_dist(A,b,lambda,ld,vec,xi,solver)

x_length = size(A,2);

%-- Parse inputs ---------------------------------------------%
if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('solver','var'); solver = []; end % if computation method not specified
%--------------------------------------------------------------%


Lpr0 = invert1d.exp_dist_lpr(ld,vec);
    % use external function to evaluate prior covariance
Lpr = lambda.*Lpr0;


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lpr],[b;sparse(x_length,1)],xi,solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end
