
function [x,D,Lx] = tikhonov(A,b,n,lambda,order,x0,solver)
% TIKHONOV  Performs inversion using various order Tikhonov regularization in two-dimensions.
% Author:   Timothy Sipkens, 2018-11-21
%
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
if ~exist('order','var') % if order not specified
    order = 1;
elseif isempty(order)
    order = 1;
end

if ~exist('x0','var') % if initial guess is not specified
    x0 = sparse(x_length,1);
elseif isempty(x0)
    x0 = sparse(x_length,1);
end

if ~exist('solver','var') % if computation method not specified
    solver = 'interior-point';
elseif isempty(solver)
    solver = 'interior-point';
end
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
switch solver
    case 'interior-point' % constrained, iterative linear least squares
        options = optimoptions('lsqlin','Algorithm','interior-point','Display','none');
        x = lsqlin([A;Lx],[b;sparse(x_length,1)],...
            [],[],[],[],x0,[],[],options);
        D = []; % not specified when using this method

    case 'trust-region-reflective'
        D = ([A;Lx]'*[A;Lx])\[A;Lx]'; % invert combined matrices to get first guess
        x0 = D*[b;Lx*zeros(x_length,1)];
        
        options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
        x = lsqlin([A;Lx],[b;sparse(x_length,1)],...
            [],[],[],[],x0,[],max(x0,0),options);

    case 'algebraic' % matrix multiplication least squares (not non-negative constrained)
        D = ([A;Lx]'*[A;Lx])\[A;Lx]'; % invert combined matrices
        x = D*[b;Lx*sparse(x_length,1)];
        
    case 'algebraic-inv' % alternate algebraic least squares (less stable than previous option)
        D = inv([A;Lx]'*[A;Lx])*[A;Lx]'; % invert combined matrices
        x = D*[b;sparse(x_length,1)];
        
end

end


function Lx = genLx1(n,x_length)
% GENLX1 Generates Tikhonov matrix for 1st order Tikhonov regularization.
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


function Lx = genLx2(n,x_length)
% GENLX2 Generates Tikhonov matrix for 2nd order Tikhonov regularization.
%-------------------------------------------------------------------------%
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


