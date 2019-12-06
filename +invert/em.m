
% EM  Performs expetation-maximization routine on data.
% Author:  Timothy Sipkens, 2019-12-06
%=========================================================================%

function [x] = em(A,b,x0)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   x0      Initial guess
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%-------------------------------------------------------------------------%


x_length = length(A(1,:));
x = x0;
x(x==0) = 1e-5; % remove initially zero values

asj = sum(A)';
asj = asj+1e-19;

for jj=1:25
    ax = A*x;
    ind_ax = ax~=0; % ignore empty rows of A

    for jj=1:x_length
        suma(jj) = sum(b(ind_ax).*A(ind_ax,jj)./ax(ind_ax));
    end
    
    x = x./asj.*suma';
end

end

