
% OPTIMIZE_EXP_DIST  Finds optimal lambda for exponential distance solver.
%=========================================================================%

function [x,lambda,out] = optimize_exp_dist(A,b,d_vec,m_vec,span,x_ex,Lex,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
%   x0      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%   Lx      Tikhonov matrix
%-------------------------------------------------------------------------%

%-- Parse inputs ---------------------------------------------------------%
if ~exist('x0','var'); x0 = []; end
if ~exist('x_ex','var'); x_ex = []; end

if ~exist('solver','var')
    x_fun = @(lambda) invert.exp_dist(A,b,d_vec,m_vec,lambda,Lex,x0);
elseif isempty(solver)
    x_fun = @(lambda) invert.exp_dist(A,b,d_vec,m_vec,lambda,Lex,x0);
else
    x_fun = @(lambda) invert.exp_dist(A,b,d_vec,m_vec,lambda,Lex,x0,solver);
end
%-------------------------------------------------------------------------%


out.lambda = logspace(log10(span(1)),log10(span(2)),70);
out.x = zeros(length(x_ex),length(out.lambda));
if ~isempty(x_ex); out.chi = zeros(length(out.lambda),1); end

disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=1:length(out.lambda)
    out.x(:,ii) = x_fun(out.lambda(ii));
    if ~isempty(x_ex); out.chi(ii) = norm(out.x(:,ii)-x_ex); end
    out.Axb(ii) = norm(A*out.x(:,ii)-b);
    
    tools.textbar(ii/length(out.lambda));
end

if ~isempty(x_ex)
    [~,ind_min] = min(out.chi);
else
    ind_min = [];
end
lambda = out.lambda(ind_min);
x = out.x(:,ind_min);

end

