function [x,iter,out] = mart_optimized(A,b,x0,iter_vec,x_ex)
% MART_OPTIMIZED Finds optimal number of iterations for MART solver.
%
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   x0      Initial guess
%   span    Span of considered MART iterations, 1x2 int
%
% Output:
%   x       MART estimate
%   iter    Optimized number of iterations
%   out     Struct containing other information
%-------------------------------------------------------------------------%

disp(' ');
disp('Optimizing MART:');
textbar(0);
disp(' ');

out.iter_vec = iter_vec;
out.x(:,1) = mart(A,b,x0,iter_vec(1));
out.chi(1) = norm(out.x(:,1)-x_ex);
textbar(1/length(iter_vec));

for ii=2:length(iter_vec)
    out.x(:,ii) = mart(A,b,out.x(:,ii-1),iter_vec(ii)-iter_vec(ii-1));
    out.chi(ii) = norm(out.x(:,ii)-x_ex);
    textbar(ii/length(iter_vec));
end

[~,ind_min] = min(out.chi);
iter = iter_vec(ind_min);
x = out.x(:,ind_min);

iter = length(iter_vec);
x = out.x(:,end);

end

