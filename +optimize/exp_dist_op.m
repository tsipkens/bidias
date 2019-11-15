
% EXP_DIST_OP  Finds optimal lambda for exponential distance solver.
%=========================================================================%

function [x,lambda,out] = exp_dist_op(A,b,d_vec,m_vec,span,x_ex,Lex,x0,solver)
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

if ~exist('Lex','var'); Lex = []; end
if isempty(Lex); Lex = speye(2); end
     % if coordinate transform is not specified

if ~exist('x0','var'); x0 = []; end % if no initial x is given
if ~exist('x_ex','var'); x_ex = []; end
%--------------------------------------------------------------%


lambda = logspace(log10(span(1)),log10(span(2)),70);

disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=length(lambda):-1:1
    out(ii).lambda = lambda(ii);
    
    [out(ii).x,~,Lpr] = invert.exp_dist(...
        A,b,d_vec,m_vec,lambda(ii),Lex,x0,solver);
    
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    out(ii).Lex = Lex;
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

if ~isempty(x_ex)
    [~,ind_min] = min([out.chi]);
else
    ind_min = [];
end
lambda = out(ind_min).lambda;
x = out(ind_min).x;

out(1).Lpr = Lpr./lambda(end); % store Lpr structure
    % to save memory, only output Lpr structure
    % Lpr for any lambda can be found using scalar multiplication
    % Gpo_inv = A'*A+Lpr'*Lpr; <- can be done is post-process

end

