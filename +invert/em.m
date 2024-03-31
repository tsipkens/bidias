
% EM  Performs expetation-maximization algorithm. 
% Assumes Poisson-distributed data.
% Author:  Timothy Sipkens, 2019-12-06
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   xi      Initial guess
%   n       Number of iterations
%
% Outputs:
%   x       Regularized estimate
%=========================================================================%

function [x] = em(A, b,xi, n)


if ~exist('n','var'); n = []; end
if isempty(n); n = 4; end

x_length = length(A(1,:));
x = xi;
x(x==0) = 1e-5; % remove initially zero values

asj = sum(A)';
asj = asj+1e-19;

for ii=1:n
    ax = A*x;
    ind_ax = ax~=0; % ignore empty rows of A

    for jj=1:x_length
        suma(jj) = sum(b(ind_ax).*A(ind_ax,jj)./ax(ind_ax));
    end
    
    x = x./asj.*suma';
end

end

