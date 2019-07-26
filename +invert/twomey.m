
% TWOMEY   Performs inversion using the iterative Twomey approach.
% Author:  Timothy Sipkens, 2018-11-21
%=========================================================================%

function [x] = twomey(A,b,x0,iter,SIGMA_fun,SIGMA,bool_textbar)
%-------------------------------------------------------------------------%
% Inputs:
%   A             Model matrix
%   b             Data
%   Lb            Cholesky factorization of inverse covariance matrix 
%   x0            Initial guess
%   iter          Max. number of iterations
%   SIGMA_fun     Function of x to be evaluated to determine convergence
%                 (Optional, default: ignore)
%   SIGMA         Value to which SIGMA_fun is compared to determine convergence
%                 (Optional, default: ignore)
%   bool_textbar  Boolean to determine whether or not to show textbar
%                 (Optional, default: 0)
%
% Outputs:
%   x             Twomey estimate
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('bool_textbar','var'); bool_textbar = []; end
if isempty(bool_textbar); bool_textbar = 0; end

if or(~exist('SIGMA','var'),~exist('SIGMA_fun','var')) % controls exit conditions
    bool_SIGMA = 0;
elseif or(isempty(SIGMA),isempty(SIGMA_fun))
    bool_SIGMA = 0;
else
    bool_SIGMA = 1;
end


%-- Start evaluation -----------------------------------------------------%
s = 1./max(A,[],2); % factor to scale data and model matrix
s_min = max(max(A,[],2))*1e-10;
s(max(A,[],2)<s_min) = 0;
A = bsxfun(@times,s,A); % scale model matrix
b = max(b.*s,0); % remove any negative data, scale data

if bool_textbar
    disp('Twomey progress:');
    tools.textbar(0); % outputs progresss
end

x = x0;
factor = 1; % factor, allows one to decrease step size in Twomey


%-- Perform Twomey iterations --------------------------------------------%
for kk=1:iter % perform multiple Twomey passes
    for ii=1:length(b) % loop through data vector for single Twomey pass
        if b(ii)~=0 % do not consider null data
            y = A(ii,:)*x;
            if y~=0 % prevents errors in dividing by zero
                X = b(ii)./y; % calculate weight factor
                C = 1+factor.*(X-1).*A(ii,:)';
                x = C.*x;
            end
        end
    end
    
    if bool_textbar; tools.textbar(kk/iter); end % outputs progresss
    
    %-- Exit conditions if SIGMA is specific (e.g. Twomey-Markowski) ----%
    if bool_SIGMA
        if SIGMA_fun(x)<SIGMA %calc_mean_sq_error(Lb*A0,x,Lb*b0)<SIGMA % average square error for cases where b~= 0
            disp(['TWOMEY: SIGMA dropped below value specified after ',num2str(kk),...
                ' iteration(s). Exiting Twomey loop.']);
            break;
        end
    end
end

end


