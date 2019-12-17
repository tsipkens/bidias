
% EXP_DIST_OPX  Finds optimal lambda for exponential distance solver.
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

function [x,lambda,out] = exp_dist_opx(A,b,d_vec,m_vec,guess,x_ex,xi,solver)

%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('xi','var'); xi = []; end % if no initial x is given
% x_ex is required
%--------------------------------------------------------------%

min_fun = @(x) norm(x-x_ex)^2;

corr_fun = @(y) min(y(4)*y(3)^2/y(2),0.98).*(y(3)^2/y(2));
Gd_fun = @(y) [(y(3)/y(2))^2,corr_fun(y);...
    corr_fun(y),y(3)^2]; % version for no correlation
    % y(2) = ratio, y(3) = ld, y(4) = corr

tic;
y0 = guess;
y1 = fminsearch(@(y) min_fun(invert.exp_dist(...
    A,b,d_vec,m_vec,y(1),Gd_fun(y),xi,solver)),...
    y0);
toc;

lambda = y1(1);
x = invert.exp_dist(...
    A,b,d_vec,m_vec,y1(1),Gd_fun(y1),xi,solver);

out.lambda = y1(1);
out.ratio = y1(2);
out.ld = y1(3);
out.corr = y1(4);
out.Gd = Gd_fun(y1);


end


