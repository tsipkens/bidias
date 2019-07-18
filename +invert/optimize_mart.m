
% OPTIMIZE_MART  Finds optimal number of iterations for MART solver.
% Author:        Timothy Sipkens, 2019-02-13
%=========================================================================%

function [x,iter,out] = optimize_mart(A,b,x0,iter_vec,x_ex)
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   x0          Initial guess
%   iter_vec    Vector fo iteration counts
%   x_ex        Exact solution
%
% Output:
%   x           MART estimate
%   iter        Optimized number of iterations
%   out         Struct containing other information
%-------------------------------------------------------------------------%

disp(' ');
disp('Optimizing MART:');
textbar(0);
disp(' ');

out.iter_vec = iter_vec;
out.x(:,1) = mart(A,b,x0,iter_vec(1));
out.chi(1) = norm(out.x(:,1)-x_ex);
tools.textbar(1/length(iter_vec));

for ii=2:length(iter_vec)
    out.x(:,ii) = mart(A,b,out.x(:,ii-1),iter_vec(ii)-iter_vec(ii-1));
    out.chi(ii) = norm(out.x(:,ii)-x_ex);
    tools.textbar(ii/length(iter_vec));
end

iter = length(iter_vec);
x = out.x(:,end);

end

