
% EXP_DIST_OPBF  Approximates optimal lambda for exponential distance solver by brute force method.
%=========================================================================%

function [x,lambda,out] = exp_dist_opbf(A,b,d_vec,m_vec,guess,x_ex,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
%   Lex     Transformation to rotate space (Optional, default is indentity matrix)
%   x0      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('x0','var'); x0 = []; end % if no initial x is given
% x_ex is required
%--------------------------------------------------------------%

min_fun = @(x) norm(x-x_ex)^2;

Gd_fun = @(y) [y(2)^2,0;0,y(3)^2]; % version for no correlation
% y(2) = sm, y(3) = sd

lambda_vec = logspace(log10(0.1*guess(1)),log10(10*guess(1)),8);
sm_vec = logspace(log10(0.1*guess(2)),log10(10*guess(2)),10);
sd_vec = logspace(log10(0.1*guess(3)),log10(10*guess(3)),10);

[vec_l,vec_m,vec_d] = ndgrid(lambda_vec,sm_vec,sd_vec);
vec_l = vec_l(:);
vec_m = vec_m(:);
vec_d = vec_d(:);

disp('Optimizing exponential distance regularization (using least-squares)...');

tools.textbar(0);
x = zeros(length(vec_l),length(x_ex));
chi = zeros(length(vec_l),1);
out(length(vec_l)).chi = [];
for ii=1:length(vec_l)
    y = [vec_l(ii),vec_m(ii),vec_d(ii)];
    x(ii,:) = invert.exp_dist(...
        A,b,d_vec,m_vec,y(1),Gd_fun(y),x0,solver);
    chi(ii) = min_fun(squeeze(x(ii,:))');
    tools.textbar(ii/length(vec_l));
    
    out.x(ii) = x(ii);
    out.chi(ii) = chi(ii);
    out(ii).lambda = vec_l(ii);
    out(ii).sm = vec_m(ii);
    out(ii).sd = vec_d(ii);
end

[~,ind_min] = min(chi);
x = out(ind_min).x;
lambda = out(ind_min).lambda;

end


