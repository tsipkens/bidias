
function [x] = twomey(A,b,x0,iter,Lb,SIGMA)
% TWOMEY Performs inversion using the iterative Twomey approach.
% Author:   Timothy Sipkens, 2018-11-21
%
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   Lb          Cholesky factorization of inverse covariance matrix 
%   x0          Initial guess
%   iter        Max. number of iterations
%
% Outputs:
%   x           Estimate
%   SIGMA_vec   Mean square error, vector over iterations
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if or(~exist('SIGMA','var'),~exist('Lb','var')) % controls exit conditions
    opt_SIGMA = 0;
elseif or(isempty(SIGMA),isempty(Lb))
    opt_SIGMA = 0;
else
    opt_SIGMA = 1;
end

b0 = b;
A0 = A;

s = 1./max(A,[],2); % factor to scale data and model matrix
s_min = max(max(A,[],2))*1e-10;
s(max(A,[],2)<s_min) = 0;
A = bsxfun(@times,s,A); % scale model matrix
b = max(b.*s,0); % remove any negative data, scale data

% fprintf('Twomey progress: ');
% textbar(0); % outputs progresss

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
    
%     y = A*x;
%     nz = and(y~=0,b~=0); % boolean of non-zero values
%     X = b(nz)./y(nz);
%     C = 1+A(nz,:)'*(X-1);
%     x = C.*x;
    
%     if mod(kk,1)==0
%     	textbar(kk/iter); % outputs progresss
%     end
    
    %-- Exit conditions if SIGMA is specific (e.g. Twomey-Markowski) ----%
    if opt_SIGMA
        if calcMeanSqErr(Lb*A0,x,Lb*b0)<SIGMA % average square error for cases where b~= 0
            disp(['TWOMEY: SIGMA dropped below value specified after ',num2str(kk),...
                ' iteration(s). Exiting Twomey loop.']);
            break;
        end
    end
end

end

