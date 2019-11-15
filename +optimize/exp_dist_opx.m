
% EXP_DIST_OPX  Finds optimal lambda for exponential distance solver.
%=========================================================================%

function [x,lambda,out] = exp_dist_opx(A,b,d_vec,m_vec,guess,x_ex,x0,solver)
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

Lex_fun = @(y2) [1/y2,0;0,2.7/y2];
tic;
disp('Optimizing exponential distance regularization (using least-squares)...');
y0 = log10(guess);
y1 = fminsearch(@(y) min_fun(invert.exp_dist(...
    A,b,d_vec,m_vec,10^y(1),Lex_fun(10^y(2)),x0,solver)),...
    y0);
toc;

lambda = y1(1);
x = invert.exp_dist(...
    A,b,d_vec,m_vec,10^y1(1),Lex_fun(10^y1(2)),x0,solver);

out.lambda = 10^y1(1);
out.s1 = 10^y1(2);
out.s2 = out.s1/2.7;
out.Lex = Lex_fun(10^y1(2));


end


% function [per] = calc_percent(ii,jj,kk,len_lam,len_s1,len_s2)
% 
% ind = (len_lam-ii+1).*(len_s1-jj+1).*(len_s2-kk+1);
% per = ind./(len_lam.*len_s1.*len_s2);
% 
% end


