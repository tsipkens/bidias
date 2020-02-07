
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
% Author:   Timothy Sipkens, 2018-11-21
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   lambda      Regularization parameter
%   order       Order of regularization -OR- 
%               pre-computed Tikhonov matrix structure
%                   (OPTIONAL, default is set by tikhonov_lpr)
%   n           The length of first dimension of solution -OR-
%               a grid, so as to support partial grids
%                   (ONLY REQUIRED when order is specified)
%   xi          Initial guess for solver
%                   (OPTIONAL, default is zeros)
%   solver      Solver (OPTIONAL, default is interior-point)
%
% Outputs:
%   x           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   Lpr0        Tikhonov matrix structure
%   Gpo_inv     Inverse of posterior covariance
%=========================================================================%

function [x,D,Lpr0,Gpo_inv] = tikhonov(A,b,lambda,order,n,xi,solver)

x_length = size(A,2);

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
    % if order not specified, use default of tikhonov_lpr

if ~exist('xi','var'); xi = []; end % if initial guess is not specified
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%


%-- Get Tikhonov smoothing matrix ----------------------------------------%
if all(size(order)==[1,1]) % if order is specified, build Lpr0
    Lpr0 = invert.tikhonov_lpr(...
        order,n,x_length);
else % is Lpr0 strucutre is provided directly
    Lpr0 = order;
end
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

