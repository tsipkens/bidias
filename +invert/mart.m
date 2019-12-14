
% MART          Vectorized inplementation of the MART algorithm
% Citation:     Gordon et al, J. Theor. Biol. 29(3), 471-481 (1970)
% Author:       Timothy Sipkens, 2018-20-12
% Adapted from: Samuel Grauer
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   iter    MART iterations     (Optional, default 10)
%   xi      Initial guess       (Optional, default ones vector)
%
% Output:
%   x       MART estimate
%=========================================================================%

function [x] = mart(A,b,xi,iter)


%-- Parse inputs -------------%
if ~exist('iter','var')
    iter = 10; % default of 10 iterations
end

if exist('xi','var')
    x = xi;
else
    x = ones(size(A,2),1); % intiate vector of ones if xi is not specified
end


%-- Get ray weights -------------%
s = 1./max(A,[],2); % factor to scale data and model matrix
s_min = max(max(A,[],2))*1e-12;
s(max(A,[],2)<s_min) = 0;
A = bsxfun(@times,s,A); % scale model matrix
b = max(b.*s,0); % remove any negative data and scale data
w = 1; % relaxation parameter, order unity


%{
%-- Note: Uncomment to switch to the normal implementation.
%-- MART iterations: Normal implementation -------------%
for kk=1:iter
    for ii = 1:length(b)
        if b(ii)~=0 % do not consider null data
            y = A(ii,:)*x; % scalar value
            if y~=0 % avoid division by zero, do not update those entries
                lam = log(b(ii)/y);
                x = x.*exp(w.*A(ii,:)'*lam);
            end
        end
    end
end
%}


%-- MART iterations: Block implementation --------------%
for kk = 1:iter
    y = A*x;
    nz = and(y~=0,b~=0); % boolean of non-zero values
    lam = log(b(nz)./y(nz));
    x = x.*exp(w.*A(nz,:)'*lam);

    if sum(isnan(x))~=0
        disp(' ');
    end
end


end
