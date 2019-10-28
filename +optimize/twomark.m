
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
if ~exist('opt_smooth','var'); opt_smooth = []; end
if isempty(opt_smooth); opt_smooth = 'Buckley'; end

if ~exist('x_ex','var'); x_ex = []; end
%-------------------------------------------------------------------------%


x_fun = @(Sf) invert.twomark(A,b,Lb,n,x0,iter,opt_smooth,Sf);

Sf = logspace(log10(span(1)),log10(span(2)),35);

disp(' ');
disp('Optimizing Twomey-Markowski smoothing:');
tools.textbar(0);
disp(' ');
for ii=length(Sf):-1:1 % loop through values of Sf
    out(ii).Sf = Sf(ii);
    out(ii).x = x_fun(out(ii).Sf);
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    disp(' ');
    disp('Optimizing Twomey-Markowski smoothing:');
    tools.textbar(0); % reinitialize textbar
    tools.textbar((length(Sf)-ii+1)/length(Sf));
    disp(' ');
end

if ~isempty(x_ex)
    [~,ind_min] = min([out.chi]);
else
    ind_min = [];
end
Sf = out(ind_min).Sf;
x = out(ind_min).x;

end


