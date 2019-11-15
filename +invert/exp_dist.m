
% EXP_DIST  Regularization based on the exponential of the inverse distance. 
% Author:   Timothy Sipkens, 2018-10-22
%=========================================================================%

function [x,D,Lpr,Gpo_inv] = exp_dist(A,b,d_vec,m_vec,lambda,Lex,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A        Model matrix
%   b        Data
%
% Outputs:
%   x        Regularized estimate
%   D        Inverse operator (x = D*[b;0])
%   Lpr      Cholesky factorization of prior covariance
%   Gpo_inv  Inverse of posterior covariance
%-------------------------------------------------------------------------%


x_length = length(A(1,:));


%-- Parse inputs ---------------------------------------------%
if ~exist('solver','var'); solver = []; end
    % if computation method not specified

if ~exist('Lex','var'); Lex = []; end
if isempty(Lex); Lex = speye(2); end
     % if coordinate transform is not specified

if ~exist('x0','var'); x0 = []; end % if no initial x is given
%--------------------------------------------------------------%


%-- Generate prior covariance matrix -----------------%
[vec_d1,vec_d2] = ndgrid(d_vec,d_vec);
[vec_m1,vec_m2] = ndgrid(m_vec,m_vec);

drm = log10(vec_m1)-log10(vec_m2);
drd = log10(vec_d1)-log10(vec_d2);
d = sqrt((drm.*Lex(1,1)+drd.*Lex(1,2)).^2+...
    (drm.*Lex(2,1)+drd.*Lex(2,2)).^2); % distance

%-- Generate prior covariance matrix --------------------------------------
Gpr = exp(-d);

Gpr_inv = inv(Gpr);
[Lpr,~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory
Lpr = lambda.*Lpr;
Lpr(abs(Lpr)<(0.0001.*max(max(abs(Lpr))))) = 0;
Lpr(d>2) = 0;
Lpr = sparse(Lpr);


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lpr],[b;sparse(x_length,1)],x0,solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end

