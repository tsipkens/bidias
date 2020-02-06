
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
% Author:   Timothy Sipkens, 2018-11-21
%-------------------------------------------------------------------------%
% Inputs:
%   A        Model matrix
%   b        Data
%   n_grid   Either (i) the length of first dimension of solution or
%            (ii) a grid, so as to support partial grids.
%   lambda   Regularization parameter
%   order    Order of regularization     (Optional, default is 1)
%   xi       Initial guess for solver    (Optional, default is zeros)
%   solver   Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x        Regularized estimate
%   D        Inverse operator (x = D*[b;0])
%   Lpr0     Tikhonov matrix
%   Gpo_inv  Inverse of posterior covariance
%=========================================================================%

function [x,D,Lpr0,Gpo_inv] = tikhonov(A,b,n_grid,lambda,order,xi,solver)

x_length = size(A,2);

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
if isempty(order); order = 1; end
    % if order not specified

if ~exist('xi','var'); xi = []; end % if initial guess is not specified
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%


%-- Generate Tikhonov smoothing matrix -----------------------------------%
Lpr0 = invert.tikhonov_lpr(x_length,n_grid,order);
Lpr = lambda.*Lpr0;


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lpr],[b;sparse(x_length,1)],xi,solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end
%=========================================================================%

