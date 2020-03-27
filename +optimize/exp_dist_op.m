
% EXP_DIST_OP  Approximates optimal lambda for exponential distance solver.
% Uses a brute force / simple discretization method over the supplied 
% 'span' for lambda.
% Author: Timothy Sipkens, 2019-12-19
% 
% Inputs:
%   A           Model matrix or Lb*A
%   b           Data or Lb*b
%   lambda      Regularization parameter
%   Gd          Mass-mobility covariance matrix, used to calculate
%               Mahalanobis distance (Optional, default: identity matrix)
%   grid_vec2	Position of elements in mobility space, grid.elements(:,2)
%   vec1        Position of elements in mass space, grid.elements(:,1)
%   x_ex        Exact distribution project to current basis (Optional, default is empty)
%   xi          Initial guess for solver    (Optional, default is ones)
%   solver      Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x           Regularized estimate
%   lambda      Semi-optimal regularization parameter
%               (against exact solution if x_ex is specified or using Bayes factor)
%   output      Output structure with information for a range of the 
%               regularization parameter
%=========================================================================%

function [x,lambda,output] = exp_dist_op(A,b,span,Gd,grid_vec2,vec1,x_ex,xi,solver,n)


%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Gd','var'); Gd = []; end
if isempty(Gd); Gd = speye(2); end
     % if coordinate transform is not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('x_ex','var'); x_ex = []; end

if ~exist('n','var'); n = []; end
if isempty(n); n = 30; end % default number of lambda entries to consider
%--------------------------------------------------------------%


lambda = logspace(log10(span(1)),log10(span(2)),n);

Lpr = invert.exp_dist_lpr(Gd,grid_vec2,vec1); % same structure throughout
[~,~,~,S1,S2] = gsvd(full(A),full(Lpr)); % pre-compute GSVD

disp('Optimizing exp. dist. regularization w.r.t. lambda...');
tools.textbar(0);
for ii=length(lambda):-1:1 % reverse loop to pre-allocate
    %-- Store case parameters ----------------------%
    output(ii).lambda = lambda(ii);
    output(ii).lm = sqrt(Gd(1,1));
    output(ii).ld = sqrt(Gd(2,2));
    output(ii).R12 = Gd(1,2)/sqrt(Gd(1,1)*Gd(2,2));
    
    %-- Perform inversion --------------------------%
    output(ii).x = invert.exp_dist(...
        A,b,lambda(ii),Gd,grid_vec2,vec1,xi,solver);
    
    %-- Store ||Ax-b|| and Euclidean error ---------%
    if ~isempty(x_ex); output(ii).eps = norm(output(ii).x-x_ex); end
    output(ii).Axb = norm(A*output(ii).x-b);
    
    %-- Compute credence, fit, and Bayes factor ----%
    [output(ii).B,output(ii).F,output(ii).C] = ...
        optimize.bayesf_precomp(A,b,output(ii).x,Lpr,output(ii).lambda,S1,S2,0);
            % bypass exp_dist_bayesf, because GSVD is pre-computed
        % optimize.exp_dist_bayesf(A,b,lambda(ii),Lpr,out(ii).x);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end


%-- Specify output ----------------------------------%
if ~isempty(x_ex)
    [~,ind_min] = min([output.eps]); % use Euclidean error
else
    [~,ind_min] = max([output.B]);
end
lambda = output(ind_min).lambda;
x = output(ind_min).x;

end

