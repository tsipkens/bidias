
% TIKHONOV_OP  Finds optimal lambda for Tikhonov solver using known distribution, x.
% Author: Timothy Sipkens, 2019-07-17
%=========================================================================%

function [x,lambda,out] = tikhonov_op(A,b,n,span,x_ex,order,xi,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
%   order   Order of regularization     (Optional, default is 1)
%   xi      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%   Lx      Tikhonov matrix
%-------------------------------------------------------------------------%

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
if ~exist('xi','var'); xi = []; end
if ~exist('x_ex','var'); x_ex = []; end
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%

lambda = logspace(log10(span(1)),log10(span(2)),70);

disp('Optimizing Tikhonov regularization:');
tools.textbar(0);
for ii=length(lambda):-1:1
    out(ii).lambda = lambda(ii);
    
    [out(ii).x,~,Lpr] = invert.tikhonov(...
        A,b,n,lambda(ii),order,xi,solver);
    
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

if ~isempty(x_ex)
    [~,ind_min] = min([out.chi]);
else
    ind_min = [];
end
lambda = out(ind_min).lambda;
x = out(ind_min).x;
out(1).ind_min = ind_min;

out(1).Lpr = Lpr; % store Lpr structure
    % to save memory, only output Lpr structure
    % Lpr for any lambda can be found using scalar multiplication
    % Gpo_inv = A'*A+Lpr'*Lpr; <- can be done is post-process

end

