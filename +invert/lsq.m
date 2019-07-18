
% LSQ      Performs least-squares or equivalent optimization.
%          For use with inversion functions like tikhonov.m 
%          and exp_dist.m
% Author:  Timothy Sipkens, 2019-07-17
%=========================================================================%

function [x,D] = lsq(A,b,solver,x0)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   x0      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
x_length = length(A(1,:));

if ~exist('x0','var'); x0 = []; end
if ~exist('solver','var'); solver = []; end

if isempty(x0); x0 = sparse(x_length,1); end % if initial guess is not specified
if isempty(solver); solver = 'interior-point'; end % if computation method not specified
%-------------------------------------------------------------------------%


%-- Perform least-squares ------------------------------------------------%
switch solver
    case 'interior-point' % constrained, iterative linear least squares
        options = optimoptions('lsqlin','Algorithm','interior-point','Display','none');
        x = lsqlin(A,b,...
            [],[],[],[],x0,[],[],options);
        D = []; % not specified when using this method

    case 'trust-region-reflective'
        D = (A'*A)\A'; % invert combined matrices to get first guess
        x0 = D*b; % Note: previously, [b;Lx*zeros(x_length,1)]
        
        options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
        x = lsqlin(A,b,...
            [],[],[],[],x0,[],max(x0,0),options);

    case 'algebraic' % matrix multiplication least squares (not non-negative constrained)
        D = (A'*A)\A'; % invert combined matrices
        x = D*b;
        
    case 'algebraic-inv' % alternate algebraic least squares (less stable than previous option)
        D = inv(A'*A)*A'; % invert combined matrices
        x = D*b;
end
%-------------------------------------------------------------------------%

end

