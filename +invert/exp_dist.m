
% EXP_DIST  Regularization based on the exponential of the distance between elements/pixels.
% Author:   Timothy Sipkens, 2018-10-22
%-------------------------------------------------------------------------%
% Inputs:
%   A        Model matrix
%   b        Data
%   grid_d   Position of elements in mobility space
%   m        Position of elements in mass space
%   lambda   Regularization parameter
%   Gd       Mass-mobility covariance matrix, used to calculate
%            Mahalanobis distance (Optional, default: identity matrix)
%   xi       Initial guess for solver (Optional, default: empty)
%   solver   Least-squares solver to use (e.g. 'interior-point', see 'lsq' function)
%
% Outputs:
%   x        Regularized estimate
%   D        Inverse operator (x = D*[b;0])
%   Lpr0     Cholesky factorization of prior covariance
%   Gpo_inv  Inverse of the posterior covariance
%=========================================================================%

function [x,D,Lpr0,Gpo_inv] = exp_dist(A,b,grid_d,m,lambda,Gd,xi,solver)


x_length = length(A(1,:));


%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Gd','var'); Gd = []; end
if isempty(Gd); Gd = speye(2); end % if not specified, use an identity matrix
if Gd(1,2)/sqrt(Gd(1,1)*Gd(2,2))>=1 % check if correlation is unphysical
    error('Correlation greater than 1.');
end

if ~exist('xi','var'); xi = []; end % if no initial x is given
%--------------------------------------------------------------%


Lpr0 = invert.exp_dist_lpr(grid_d,m,Gd);
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
