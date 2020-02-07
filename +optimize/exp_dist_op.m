
% EXP_DIST_OP  Finds optimal lambda for exponential distance solver.
% Author: Timothy Sipkens, 2019-12-19
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
%   Lex     Transformation to rotate space (Optional, default is indentity matrix)
%   xi      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%=========================================================================%

function [x,lambda,out] = exp_dist_op(A,b,span,Gd,d_vec,m_vec,x_ex,xi,solver)


%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Gd','var'); Gd = []; end
if isempty(Gd); Gd = speye(2); end
     % if coordinate transform is not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('x_ex','var'); x_ex = []; end
%--------------------------------------------------------------%


lambda = logspace(log10(span(1)),log10(span(2)),30);

disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=length(lambda):-1:1
    %-- Store case parameters ----------------------%
    out(ii).lambda = lambda(ii);
    out(ii).lm = sqrt(Gd(1,1));
    out(ii).ld = sqrt(Gd(2,2));
    out(ii).R12 = Gd(1,2)/sqrt(Gd(1,1)*Gd(2,2));
    
    %-- Perform inversion --------------------------%
    [out(ii).x,~,Lpr] = invert.exp_dist(...
        A,b,lambda(ii),Gd,d_vec,m_vec,xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error ---------%
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor ----%
    [out(ii).B,out(ii).F,out(ii).C] = ...
        optimize.exp_dist_bayesf(A,b,lambda(ii),Lpr,out(ii).x);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

if ~isempty(x_ex)
    [~,ind_min] = min([out.chi]);
else
    ind_min = [];
end
lambda = out(ind_min).lambda;
x = out(ind_min).x;

end

