
% EXP_DIST  Regularization based on the exponential of the inverse distance. 
% Author:   Timothy Sipkens, 2018-10-22
%=========================================================================%

function [x,D,Lx] = exp_dist(A,b,d_vec,m_vec,lambda,Lex,x0,solver,sigma)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%   Lx      Cholesky factorization of posterior covariance
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

Gx = exp(-d);
if exist('sigma','var') % incorporate structure into covariance, if specified
    for ii=1:x_length
        for jj=1:x_length
            Gx(ii,jj) = Gx(ii,jj).*sigma(ii).*sigma(jj);
        end
    end
end
Gx(Gx<(0.05.*mean(mean(Gx)))) = 0; % remove any entries below thershold
Gx = Gx./max(max(Gx)); % normalize matrix structure

Gxi = inv(Gx);
Lx = chol(Gxi);
Lx = lambda.*Lx./max(max(Lx));
Lx(abs(Lx)<(0.01.*mean(mean(abs(Lx))))) = 0;
Lx = sparse(Lx);


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lx],[b;sparse(x_length,1)],solver,x0);

end

