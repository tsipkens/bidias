
% TWOMEY_OP  Finds optimal number of iterations for MART solver.
%
% Inputs:
%   A           Model matrix
%   b           Data
%   xi          Initial guess
%   iter_vec    Span of considered Twomey iterations, 1x2 int
%
% Output:
%   x           MART estimate
%   iter        Optimized number of iterations
%   output      Output structure containing other information
%=========================================================================%

function [x,iter,out] = twomey_op(A,b,xi,iter_vec,x_ex)

disp(' ');
disp('Optimizing Twomey:');
tools.textbar(0);
disp(' ');

out(length(iter_vec)).iter_vec = [];
    % initialize size of output structure

out(1).iter_vec = iter_vec(1);
out(1).x = invert.twomey(A,b,xi,iter_vec(1));
out(1).chi = norm(out(1).x-x_ex);
tools.textbar(1/length(iter_vec));

for ii=2:length(iter_vec)
    out(ii).iter_vec = ii;
    out(ii).x = invert.twomey(A,b,out(ii-1).x,1);
    out(ii).chi = norm(out(ii).x-x_ex);
    tools.textbar(ii/length(iter_vec));
end

iter = length(iter_vec);
x = out(end).x;

end

