
% TIKHONOV_OP  Finds optimal lambda for Tikhonov solver using known distribution, x.
% Author: Timothy Sipkens, 2019-07-17
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
%=========================================================================%

function [x,lambda,out] = tikhonov_op(A,b,n,span,x_ex,order,xi,solver)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
if ~exist('xi','var'); xi = []; end
if ~exist('x_ex','var'); x_ex = []; end
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%

lambda = logspace(log10(span(1)),log10(span(2)),70);

[~,~,Lpr] = invert.tikhonov(...
    A,b,n,lambda(end),order,xi,solver); % just used to get Lpr
[~,~,~,S1,S2] = gsvd(full(A),full(Lpr));
    % pre-compute gsvd for Bayes factor calculation

disp('Optimizing Tikhonov regularization:');
tools.textbar(0);
for ii=length(lambda):-1:1
    out(ii).lambda = lambda(ii); % store regularization parameter
    
    %-- Perform inversion --%
    [out(ii).x,~,Lpr] = invert.tikhonov(...
        A,b,n,lambda(ii),order,xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [out(ii).B,out(ii).F,out(ii).C] = ...
        optimize_b.tikhonov_bayesf(A,b,out(ii).x,Lpr,lambda(ii),order,...
        S1,S2);
    
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

