
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
% Author:   Timothy Sipkens, 2018-11-21
%=========================================================================%

function [x,D,Lx] = tikhonov(A,b,n,lambda,order,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   order   Order of regularization     (Optional, default is 1)
%   x0      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%   Lx      Tikhonov matrix
%-------------------------------------------------------------------------%

x_length = length(A(1,:));

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
if ~exist('x0','var'); x0 = []; end % if initial guess is not specified
if ~exist('solver','var'); solver = []; end

if isempty(order); order = 1; end % if order not specified
% if isempty(x0); x0 = sparse(x_length,1); end % if initial guess is not specified
if isempty(solver); solver = 'interior-point'; end % if computation method not specified
%-------------------------------------------------------------------------%


%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
        Lx = -lambda.*speye(x_length);
    case 1 % 1st order Tikhonov
        Lx = lambda.*genLx1(n,x_length);
    case 2 % 2nd order Tikhonov
        Lx = lambda.*genLx2(n,x_length);
    case 10
        Lx = -0.1.*lambda.*speye(x_length)+lambda.*genLx1(n,x_length);
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lx],[b;sparse(x_length,1)],solver,x0);

end
%=========================================================================%


%== GENLX1 ===============================================================%
%   Generates Tikhonov matrix for 1st order Tikhonov regularization.
function Lx = genLx1(n,x_length)
%-------------------------------------------------------------------------%
% Inputs:
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   Lx      Tikhonov matrix
%-------------------------------------------------------------------------%

% Dx = speye(n);
% Dx = spdiag(-ones(n,1),1,Dx);
% Dx = kron(Dx,speye(x_length/n));

Lx = -eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        Lx(jj,jj+1) = 0.5;
    else % if on right edge
        Lx(jj,jj) = Lx(jj,jj)+0.5;
    end
    
    if jj<=(x_length-n)
        Lx(jj,jj+n) = 0.5;
    else % if on bottom
        Lx(jj,jj) = Lx(jj,jj)+0.5;
    end
end
Lx = sparse(Lx);
 
end
%=========================================================================%


%== GENLX2 ===============================================================%
%   Generates Tikhonov matrix for 2nd order Tikhonov regularization.
function Lx = genLx2(n,x_length)
% Inputs:
%   n           Length of first dimension of solution
%   x_length    Length of x vector
%
% Outputs:
%   Lx      Tikhonov matrix
%-------------------------------------------------------------------------%

Lx = -eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        Lx(jj,jj+1) = 0.25;
    else
        Lx(jj,jj) = Lx(jj,jj)+0.25;
    end
    
    if ~(mod(jj-1,n)==0)
        Lx(jj,jj-1) = 0.25;
    else
        Lx(jj,jj) = Lx(jj,jj)+0.25;
    end
    
    if jj>n
        Lx(jj,jj-n) = 0.25;
    else
        Lx(jj,jj) = Lx(jj,jj)+0.25;
    end
    
    if jj<=(x_length-n)
        Lx(jj,jj+n) = 0.25;
    else
        Lx(jj,jj) = Lx(jj,jj)+0.25;
    end
end
Lx = sparse(Lx);

end
%=========================================================================%

