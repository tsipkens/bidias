
% OPTIMIZE_TWOMARK  Finds optimal smoothing for Twomey-Markowski solver using known distribution, x.
% Author:           Timothy Sipkens, 2018-12-20
%=========================================================================%

function [x,Sf,out] = optimize_twomark(A,b,Lb,n,x0,iter,span,x_ex,opt_smooth)
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   Lb          Cholesky factorization of inverse covariance matrix
%   n           Length of first dimension of solution, used in smoothing
%   x0          Initial guess
%   iter        Max. number of iterations of Twomey_Markowski algorithm
%   span        Range for 1/Sf, two entry vector
%   x_ex        Exact distribution project to current basis
%   opt_smooth  Type of smoothing to apply      (Optional, default is 'Buckley')
%
% Outputs:
%   x           Estimate
%   Sf          Estimate of optimized value of Sf
%   out         Struct containing detailed information about optimization
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('opt_smooth','var')
    opt_smooth = 'Buckley';
elseif isempty(opt_smooth)
    opt_smooth = 'Buckley';
end
%-------------------------------------------------------------------------%


x_fun = @(Sf) twomey_markowski(A,b,Lb,n,x0,iter,opt_smooth,Sf);

out.Sf = 1./logspace(log10(span(1)),log10(span(2)),25);
out.x = zeros(length(x_ex),length(out.Sf));
out.chi = zeros(length(out.Sf),1);

disp(' ');
disp('Optimizing Twomey-Markowski smoothing:');
textbar(0);
disp(' ');
for ii=1:length(out.Sf)
    out.x(:,ii) = x_fun(out.Sf(ii));
    out.chi(ii) = norm(out.x(:,ii)-x_ex);
    out.Axb(ii) = norm(A*out.x(:,ii)-b);
    
    disp(' ');
    disp('Optimizing Twomey-Markowski smoothing:');
    textbar(ii/length(out.Sf));
    disp(' ');
end

[~,ind_min] = min(out.chi);
Sf = out.Sf(ind_min);
x = out.x(:,ind_min);

end


