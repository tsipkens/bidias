
function [x] = mart(A,b,x0,iter)
% MART Fast inplementation of the MART algorithm
% Citation:     Gordon et al, J. Theor. Biol. 29(3), 471-481 (1970)
% Author:       Samuel Grauer, 2018-10-01
% Modified by:  Timothy Sipkens, 2018-20-12
%
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   iter    MART iterations     (Optional, default 10)
%   x0      Initial guess       (Optional, default ones vector)
%
% Output:
%   x       MART estimate
%-------------------------------------------------------------------------%


%-- Parse inputs -------------%
if ~exist('iter','var')
    iter = 10; % default of 10 iterations
end

if exist('x0','var')
    x = x0;
else
    x = ones(size(A,2),1); % intiate vector of ones if x0 is not specified
end


%-- Get ray weights -------------%
% s = 1./sum(A)';
s = 1./max(A,[],2); % factor to scale data and model matrix
s_min = max(max(A,[],2))*1e-12;
s(max(A,[],2)<s_min) = 0;
A = bsxfun(@times,s,A); % scale model matrix
b = max(b.*s,0); % remove any negative data and scale data
w = 1; % relaxation parameter, order unity


%-- MART iterations: Normal implementation -------------%
% for kk=1:iter
%     for ii = 1:length(b)
%         if b(ii)~=0 % do not consider null data
%             y = A(ii,:)*x; % scalar value
%             if y~=0 % avoid division by zero, do not update those entries
%                 lam = log(b(ii)/y);
%                 x = x.*exp(w.*A(ii,:)'*lam);
%             end
%         end
%     end
% end


%-- MART iterations: Block implementation --------------%
for kk = 1:iter
    y = A*x;
    nz = and(y~=0,b~=0); % boolean of non-zero values
    lam = log(b(nz)./y(nz));
    x = x.*exp(w.*A(nz,:)'*lam);
end


%-- Apply non-negativity constraint -------------%
x = max(0,real(x));

end


