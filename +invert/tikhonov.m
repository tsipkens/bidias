
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
% Author:   Timothy Sipkens, 2018-11-21
%=========================================================================%

function [x,D,Lpr0,Gpo_inv] = tikhonov(A,b,n,lambda,order,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A        Model matrix
%   b        Data
%   n        Length of first dimension of solution
%   lambda   Regularization parameter
%   order    Order of regularization     (Optional, default is 1)
%   x0       Initial guess for solver    (Optional, default is zeros)
%   solver   Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x        Regularized estimate
%   D        Inverse operator (x = D*[b;0])
%   Lpr0     Tikhonov matrix
%   Gpo_inv  Inverse of posterior covariance
%-------------------------------------------------------------------------%

x_length = length(A(1,:));

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
if isempty(order); order = 1; end
    % if order not specified

if ~exist('x0','var'); x0 = []; end % if initial guess is not specified
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%


%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
        Lpr0 = -speye(x_length);
    case 1 % 1st order Tikhonov
        Lpr0 = genL1(n,x_length);
    case 2 % 2nd order Tikhonov
        Lpr0 = genL2(n,x_length);
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end
Lpr = lambda.*Lpr0;


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lpr],[b;sparse(x_length,1)],x0,solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end
%=========================================================================%


%== GENL1 ================================================================%
%   Generates Tikhonov matrix for 1st order Tikhonov regularization.
function L = genL1(n,x_length)
%-------------------------------------------------------------------------%
% Inputs:
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   L       Tikhonov matrix
%-------------------------------------------------------------------------%

% Dx = speye(n);
% Dx = spdiag(-ones(n,1),1,Dx);
% Dx = kron(Dx,speye(x_length/n));

L = -eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        L(jj,jj+1) = 0.5;
    else % if on right edge
        L(jj,jj) = L(jj,jj)+0.5;
    end
    
    if jj<=(x_length-n)
        L(jj,jj+n) = 0.5;
    else % if on bottom
        L(jj,jj) = L(jj,jj)+0.5;
    end
end
L = sparse(L);
 
end
%=========================================================================%


%== GENL2 ================================================================%
%   Generates Tikhonov matrix for 2nd order Tikhonov regularization.
function L = genL2(n,x_length)
% Inputs:
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   L       Tikhonov matrix
%-------------------------------------------------------------------------%

L = -eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        L(jj,jj+1) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end
    
    if ~(mod(jj-1,n)==0)
        L(jj,jj-1) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end
    
    if jj>n
        L(jj,jj-n) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end
    
    if jj<=(x_length-n)
        L(jj,jj+n) = 0.25;
    else
        L(jj,jj) = L(jj,jj)+0.25;
    end
end
L = sparse(L);

end
%=========================================================================%

