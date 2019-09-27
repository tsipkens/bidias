
% EXP_DIST  Regularization based on the exponential of the inverse distance. 
% Author:   Timothy Sipkens, 2018-10-22
%=========================================================================%

function [x,D,Lpr,Gpo_inv] = exp_dist(A,b,d_vec,m_vec,lambda,Lex,x0,solver,sigma)
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
if isempty(solver); solver = 'interior-point'; end
    % if computation method not specified

if ~exist('Lex','var'); Lex = []; end
if isempty(solver); Lex = speye(2); end
     % if coordinate transform is not specified

if ~exist('x0','var'); x0 = []; end % if no initial x is given
%--------------------------------------------------------------%


%-- Generate prior covariance matrix -----------------%
[vec_d1,vec_d2] = ndgrid(d_vec,d_vec);
[vec_m1,vec_m2] = ndgrid(m_vec,m_vec);

d1 = log(vec_m1)-log(vec_m2);
d2 = log(vec_d1)-log(vec_d2);
d = sqrt((d1.*Lex(1,1)+d2.*Lex(1,2)).^2+(d1.*Lex(2,1)+d2.*Lex(2,2)).^2); % distance

%-- Generate prior covariance matrix --------------------------------------
Gpr = exp(-d);
if exist('sigma','var') % incorporate structure into covariance, if specified
    for ii=1:x_length
        for jj=1:x_length
            Gpr(ii,jj) = Gpr(ii,jj).*sigma(ii).*sigma(jj);
        end
    end
end
Gpr(Gpr<(0.05.*mean(mean(Gpr)))) = 0; % remove any entries below thershold
Gpr = Gpr./max(max(Gpr)); % normalize matrix structure

Gpr_inv = inv(Gpr);
Lpr = chol(Gpr_inv);
clear Gpr_inv; % to save memory
Lpr = lambda.*Lpr./max(max(Lpr));
Lpr(abs(Lpr)<(0.01.*mean(mean(abs(Lpr))))) = 0;
Lpr = sparse(Lpr);


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lpr],[b;sparse(x_length,1)],solver,x0);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end

