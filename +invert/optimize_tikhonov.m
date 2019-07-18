
% OPTIMIZE_TIKHONOV  Finds optimal lambda for Tikhonov solver using known distribution, x.
%=========================================================================%

function [x,lambda,out] = optimize_tikhonov(A,b,n,span,x_ex,order,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
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
if isempty(order); order = 1; end

if ~exist('x0','var'); x0 = []; end

if ~exist('solver','var')
    x_fun = @(lambda) tikhonov(A,b,n,lambda,order,x0);
elseif isempty(solver)
    x_fun = @(lambda) tikhonov(A,b,n,lambda,order,x0);
else
    x_fun = @(lambda) tikhonov(A,b,n,lambda,order,x0,solver);
end
%-------------------------------------------------------------------------%


out.lambda = logspace(log10(span(1)),log10(span(2)),70);
out.x = zeros(length(x_ex),length(out.lambda));
out.chi = zeros(length(out.lambda),1);

disp('Optimizing Tikhonov regularization:');
tools.textbar(0);
for ii=1:length(out.lambda)
    out.x(:,ii) = x_fun(out.lambda(ii));
    out.chi(ii) = norm(out.x(:,ii)-x_ex);
    out.Axb(ii) = norm(A*out.x(:,ii)-b);
    
    tools.textbar(ii/length(out.lambda));
end

[~,ind_min] = min(out.chi);
lambda = out.lambda(ind_min);
x = out.x(:,ind_min);

end

