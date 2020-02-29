
% TIKHONOV_OP  Finds optimal lambda for Tikhonov solver using known distribution, x.
% Author: Timothy Sipkens, 2019-07-17
% 
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
%   lambda  Semi-optimal regularization parameter
%           (against exact solution if x_ex is specified or using Bayes factor)
%   output  Output structure with information for a range of the regularization parameter
%=========================================================================%

function [x,lambda,output] = tikhonov_op(A,b,span,order,n_grid,x_ex,xi,solver,n)

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
for ii=length(lambda):-1:1 % reverse loop to pre-allocate
    output(ii).lambda = lambda(ii); % store regularization parameter
    
    %-- Perform inversion --%
    [output(ii).x,~,Lpr0] = invert.tikhonov(...
        A,b,lambda(ii),Lpr0,[],xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    if ~isempty(x_ex); output(ii).chi = norm(output(ii).x-x_ex); end
    output(ii).Axb = norm(A*output(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [output(ii).B,output(ii).F,output(ii).C] = ...
        optimize.bayesf_precomp(...
        A,b,output(ii).x,Lpr0,lambda(ii),S1,S2,order);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

if ~isempty(x_ex) % if exact solution is supplied
    [~,ind_min] = min([output.chi]);
else
    [~,ind_min] = max([output.B]);
end
lambda = output(ind_min).lambda;
x = output(ind_min).x;
output(1).ind_min = ind_min;

output(1).Lpr = Lpr0; % store Lpr structure
    % to save memory, only output Lpr structure
    % Lpr for any lambda can be found using scalar multiplication
    % Gpo_inv = A'*A+Lpr'*Lpr; <- can be done is post-process

end

