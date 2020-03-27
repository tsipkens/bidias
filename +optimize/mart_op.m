
% MART_OP  Finds optimal number of iterations for MART solver.
% Author:  Timothy Sipkens, 2019-02-13
% 
% Inputs:
%   A           Model matrix
%   b           Data
%   xi          Initial guess
%   iter_vec    Vector fo iteration counts
%   x_ex        Exact solution
%
% Output:
%   x           MART estimate
%   iter        Optimized number of iterations
%   output      Output structure containing other information
%=========================================================================%

function [x,iter,output] = mart_op(A,b,xi,iter_vec,x_ex)

disp(' ');
disp('Optimizing MART:');
tools.textbar(0);
disp(' ');

output(length(iter_vec)).iter_vec = [];
    % initialize size of output structure

output(1).iter_vec = iter_vec(1);
output(1).x = invert.mart(A,b,xi,iter_vec(1));
output(1).eps = norm(output(1).x-x_ex); % Euclidean error
tools.textbar(1/length(iter_vec));

for ii=2:length(iter_vec)
    output(ii).iter_vec = ii;
    output(ii).x = invert.mart(A,b,output(ii-1).x,iter_vec(ii)-iter_vec(ii-1));
    output(ii).eps = norm(output(ii).x-x_ex); % Euclidean error
    
    if any(isnan(output(ii).x)); break; end
        % if NaN is encountered exit function (warning thrown my MART)
    
    tools.textbar(ii/length(iter_vec));
end

iter = ii;
x = output(end).x;

end

