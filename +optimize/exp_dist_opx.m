
% EXP_DIST_OPX  Uses fminsearch to find optimal regularization parameters
% 
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
%   lambda  Optimized regularization parameter
%   output  Output structure containing other information
%=========================================================================%

function [x,lambda,output] = exp_dist_opx(A,b,d_vec,m_vec,guess,x_ex,xi,solver)

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
    A,b,y(1),Gd_fun(y),d_vec,m_vec,xi,solver)),...
    y0);
toc;

lambda = y1(1);
x = invert.exp_dist(...
    A,b,y1(1),Gd_fun(y1),d_vec,m_vec,xi,solver);

output.lambda = y1(1);
output.ratio = y1(2);
output.ld = y1(3);
output.corr = y1(4);
output.Gd = Gd_fun(y1);


end


