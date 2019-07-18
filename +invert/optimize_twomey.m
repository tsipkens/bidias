
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
textbar(0);
disp(' ');

out.iter_vec = iter_vec;
out.x(:,1) = mart(A,b,x0,iter_vec(1));
out.chi(1) = norm(out.x(:,1)-x_ex);
tools.textbar(1/length(iter_vec));

for ii=2:length(iter_vec)
    out.x(:,ii) = invert.twomey(A,b,out.x(:,ii-1),1);
    out.chi(ii) = norm(out.x(:,ii)-x_ex);
    tools.textbar(ii/length(iter_vec));
end

[~,ind_min] = min(out.chi);
iter = iter_vec(ind_min);
x = out.x(:,ind_min);

iter = length(iter_vec);
x = out.x(:,end);

end

