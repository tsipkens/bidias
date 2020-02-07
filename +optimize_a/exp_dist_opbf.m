
% EXP_DIST_OPBF  Approximates optimal lambda for exponential distance solver by brute force method.
% Author: Timothy Sipkens, 2019-12
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

function [x,lambda,out] = exp_dist_opbf(A,b,d_vec,m_vec,x_ex,xi,solver)


%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
% x_ex is required
%--------------------------------------------------------------%

Gd_fun = @(y) [(y(3)/y(2))^2,y(4)*y(3)^2/y(2);y(4)*y(3)^2/y(2),y(3)^2];
    % y(2) = ratio, y(3) = ld, y(4) = corr

lambda = logspace(log10(1),log10(2),5);
ratio = logspace(log10(1/4),log10(1/2),5); % ratio = ld/lm
ld = logspace(log10(log10(1.5)),log10(log10(2.2)),5);
corr = linspace(0.89,0.99,5);

[vec_lambda,vec_ratio,vec_ld,vec_corr] = ...
    ndgrid(lambda,ratio,ld,corr);
vec_lambda = vec_lambda(:);
vec_ratio = vec_ratio(:);
vec_ld = vec_ld(:);
vec_corr = vec_corr(:);

tools.textbar(0);
out(length(vec_lambda)).chi = [];
for ii=1:length(vec_lambda)
    y = [vec_lambda(ii),vec_ratio(ii),vec_ld(ii),vec_corr(ii)];
    
    out(ii).x = invert.exp_dist(...
        A,b,y(1),Gd_fun(y),d_vec,m_vec,xi,solver);
    out(ii).chi = norm(out(ii).x-x_ex);
    
    out(ii).lambda = vec_lambda(ii);
    out(ii).ratio = vec_ratio(ii);
    out(ii).ld = vec_ld(ii);
    out(ii).corr = vec_corr(ii);
    
    [out(ii).B,out(ii).F,out(ii).C] = ...
        optimize_b.exp_dist_bayesf(A,b,out(ii).lambda,Lpr,out(ii).x);
    
    tools.textbar(ii/length(vec_lambda));
end

tools.textbar(1);

[~,ind_min] = min([out(ii).chi]);
x = [];%out(ind_min).x;
lambda = [];%out(ind_min).lambda;

end


