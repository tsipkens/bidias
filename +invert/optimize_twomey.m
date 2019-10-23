
% OPTIMIZE_TWOMEY  Finds optimal number of iterations for MART solver.
%=========================================================================%

function [x,iter,out] = optimize_twomey(A,b,x0,iter_vec,x_ex)
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   x0          Initial guess
%   iter_vec    Span of considered Twomey iterations, 1x2 int
%
% Output:
%   x           MART estimate
%   iter        Optimized number of iterations
%   out         Struct containing other information
%-------------------------------------------------------------------------%

disp(' ');
disp('Optimizing Twomey:');
tools.textbar(0);
disp(' ');

out(length(iter_vec)).chi = [];
    % initialize size of output structure

out(1).iter_vec = iter_vec(1);
out(1).x = invert.twomey(A,b,x0,iter_vec(1));
out(1).chi = norm(out(1).x-x_ex);
tools.textbar(1/length(iter_vec));

for ii=2:length(iter_vec)
    out(ii).x = invert.twomey(A,b,out(ii).x,1);
    out(ii).chi = norm(out(ii).x-x_ex);
    tools.textbar(ii/length(iter_vec));
end

if ~isempty(x_ex)
    [~,out.ind_min] = min([out.chi]);
else
    out.ind_min = [];
end
out.x_opt = out(out.ind_min).x; % store optimal solution in out struct


%-- Ouput solution at last step (for compatibility) -------%
iter = length(iter_vec);
x = out(end).x;

end

