
% TIKHONOV_OP  Finds optimal lambda for Tikhonov solver using known distribution, x.
% Author: Timothy Sipkens, 2019-07-17
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n_grid  Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
%   order   Order of regularization     (Optional, default is 1)
%   xi      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%   n       Number of lambda entries    (Optional, default is 70)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%   Lx      Tikhonov matrix
%=========================================================================%

function [x,lambda,out] = tikhonov_op(A,b,span,order,n_grid,x_ex,xi,solver,n)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end
if ~exist('xi','var'); xi = []; end
if ~exist('x_ex','var'); x_ex = []; end
if ~exist('solver','var'); solver = []; end

if ~exist('n','var'); n = []; end
if isempty(n); n = 70; end % default number of lambda entries to consider
%-------------------------------------------------------------------------%

lambda = logspace(log10(span(1)),log10(span(2)),n);
x_length = size(A,2);

Lpr0 = invert.tikhonov_lpr(order,n_grid,x_length); % get Tikhonov matrix

disp('Pre-computing GSV...');
[~,~,~,S1,S2] = gsvd(full(A),full(Lpr0));
    % pre-compute gsvd for Bayes factor calculation
disp('Complete.');
disp(' ');

disp('Optimizing Tikhonov regularization:');
tools.textbar(0);
for ii=length(lambda):-1:1
    out(ii).lambda = lambda(ii); % store regularization parameter
    
    %-- Perform inversion --%
    [out(ii).x,~,Lpr0] = invert.tikhonov(...
        A,b,lambda(ii),Lpr0,[],xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [out(ii).B,out(ii).F,out(ii).C] = ...
        optimize.bayesf_precomp(...
        A,b,out(ii).x,Lpr0,lambda(ii),S1,S2,order);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

if ~isempty(x_ex) % if exact solution is supplied
    [~,ind_min] = min([out.chi]);
else
    ind_min = [];
end
lambda = out(ind_min).lambda;
x = out(ind_min).x;
out(1).ind_min = ind_min;

out(1).Lpr = Lpr0; % store Lpr structure
    % to save memory, only output Lpr structure
    % Lpr for any lambda can be found using scalar multiplication
    % Gpo_inv = A'*A+Lpr'*Lpr; <- can be done is post-process

end

