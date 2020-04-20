
% TIKHONOV_OP  Finds optimal lambda for Tikhonov solver using known distribution, x.
% Author: Timothy Sipkens, 2019-07-17
% 
% Inputs:
%   A       Model matrix
%   b       Data
%   n_grid  Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   order   Order of regularization     (Optional, default is 1)
%   n       Number of lambda entries    (Optional, default is 70)
%
% Outputs:
%   x       Regularized estimate
%   lambda  Semi-optimal regularization parameter
%           (against exact solution if x_ex is specified or using Bayes factor)
%   output  Output structure with information for a range of the regularization parameter
%   Gpo_inv     Inverse of posterior covariance
%=========================================================================%

function [x,lambda,output,Gpo_inv] = tikhonov_op(A,b,span,order,n)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order','var'); order = []; end

if ~exist('n','var'); n = []; end
if isempty(n); n = 70; end % default number of lambda entries to consider
%-------------------------------------------------------------------------%

lambda = logspace(log10(span(1)),log10(span(2)),n);
x_length = size(A,2);

Lpr0 = invert1d.tikhonov_lpr(order,x_length); % get Tikhonov matrix

disp('Pre-computing generalized SVD...');
[~,~,~,S1,S2] = gsvd(full(A),full(Lpr0));
    % pre-compute gsvd for Bayes factor calculation
disp('Complete.');
disp(' ');

disp('Optimizing Tikhonov regularization w.r.t lambda...');
tools.textbar(0);
for ii=length(lambda):-1:1 % reverse loop to pre-allocate
    output(ii).lambda = lambda(ii); % store regularization parameter
    
    %-- Perform inversion --%
    [output(ii).x,~,Lpr0] = invert.tikhonov(...
        A,b,lambda(ii),Lpr0);
    
    %-- Store ||Ax-b|| and Euclidean error --%
    output(ii).Axb = norm(A*output(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor --%
    [output(ii).B,output(ii).F,output(ii).C] = ...
        optimize.bayesf_precomp(...
        A,b,output(ii).x,Lpr0,lambda(ii),S1,S2,order);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

[~,ind_min] = max([output.B]); % use Bayes factor

lambda = output(ind_min).lambda;
x = output(ind_min).x;
output(1).ind_min = ind_min;

output(1).Lpr = Lpr0; % store Lpr structure
    % to save memory, only output Lpr structure
    % Lpr for any lambda can be found using scalar multiplication
    % Gpo_inv = A'*A+Lpr'*Lpr; <- can be done is post-process

if nargout>=4
    Gpo_inv = A'*A+lambda^2.*(Lpr0'*Lpr0);
end

end

