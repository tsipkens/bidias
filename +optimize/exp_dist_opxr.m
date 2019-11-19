
% EXP_DIST_OPXR  Finds optimal lambda for exponential distance, rotated solver.
%=========================================================================%

function [x,lambda,out] = exp_dist_opxr(A,b,d_vec,m_vec,guess,x_ex,x0,solver)
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

thresh = @(y) max(y(2),y(4)*y(3)-1e-5); % to prevent correlation >= 1
Gd_fun = @(y) [thresh(y)^2,y(3)^2*y(4);...
    y(3)^2*y(4),y(3)^2];
% y(2) = sm, y(3) = sd, y(4) = Dm

tic;
disp('Optimizing exponential distance regularization (using least-squares)...');
y0 = guess;
y1 = fminsearch(@(y) min_fun(invert.exp_dist(...
    A,b,d_vec,m_vec,y(1),Gd_fun(y),x0,solver)),...
    y0);
toc;

y1(2) = thresh(y1);
lambda = y1(1);
x = invert.exp_dist(...
    A,b,d_vec,m_vec,10^y1(1),Gd_fun(y1),x0,solver);

out.lambda = y1(1);
out.s1 = y1(2);
out.s2 = y1(3);
out.Dm = y1(4);
out.Lex = Gd_fun(y1);


end


