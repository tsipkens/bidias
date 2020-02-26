
% GET_NOISE Simulates Poisson-Gaussian noise in signals.
% Author: Timothy Sipkens, 2019-12-08
% 
% Modeled after Sipkens et al., Appl. Opt. (2017).
% DOI: https://doi.org/10.1364/AO.56.008436
%-------------------------------------------------------------------------%
% Inputs:
%   b0      Uncorrupted data vector
% 
%   n_tot   Factor by which the data was normalized
%           Note that get_noise(n_tot.*b0)./n_tot; and
%           add_noise(b0,n_tot); should give identical results. 
%               (Optional, dafault: n_tot = 1)
% 
%   gam0    Percentage of the maximum signal used for the background
%           Gaussian noise. 
%               (Optional, default: gam0 = 1e-4, i.e. 0.01% of the peak)
% 
%   f_apx   Flag for whether to use a Gaussian approximation for Poisson
%           noise. Note that Lb is the Cholesky factorization of the
%           inverse of the covariance matrix for the Gaussian
%           approximation. Though, this covariance will 
%           well-approximates the covariance information for the
%           true Poisson case. 
%               (Optional, default: f_apx = 1;)
%=========================================================================%

function [b,Lb] = get_noise(b0,n_tot,gam0,f_apx)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('n_tot','var'); n_tot = []; end
if isempty(n_tot); n_tot = 1; end

if ~exist('gam0','var'); gam0 = []; end
if isempty(gam0); gam0 = 1e-4; end

if ~exist('f_apx','var'); f_apx = []; end
if isempty(f_apx); f_apx = 1; end
%-------------------------------------------------------------------------%


n_b = length(b0); % length of the data vector



%== POISSON NOISE ========================================================%
%   Poisson noise is approximated as Gaussian, with a standard
%   deviation of sqrt of the number of counts. If the data is unscaled
%   norm_fact = 1. If the data is scaled, which generally improves
%   inference, theta accounts for the normalization factor. As norm_fact
%   increases, the noise level drops, a consequence of their being more
%   counts in the data than the scaled b0 indicates.
theta = 1/n_tot;
sig_pois = sqrt(theta.*b0);
% equivalent to: sig_pois = sqrt(n_tot.*b0)./n_tot;
%     such that: b/n_tot = b0/n_tot + sig_pois0/n_tot



%== ADDITIVE GAUSSIAN NOISE ==============================================%
%   Take additive Gaussian noise as (per_gaus*100)% of peak signal.
%   This represents the minimum noise level and is to model
%   background source, such as electronic noise in the CPC
%   and fluctutions in the aerosol.
%   As a second note, increasing this quantity tends
%   to result in poor reconstruction in background but better
%   reconstruction of the peak of the distribution.
gamma = max(sig_pois)*gam0;
sig_gaus = gamma;



%== CALCULATE COVARIANCE INFORMATION =====================================%
%   Standard deviation of the combined noise.
%   The standard deviation of the sum of two idependent, normal random
%   variables is the sqrt of the sum of their squares.
sig = sqrt(sig_pois.^2+sig_gaus^2);

Lb = sparse(1:n_b,1:n_b,1./sig,n_b,n_b);
    % Cholesky factorization of the inverse covariance matrix
    % As the noise is independnet, this is a diagonal matrix with
    % 1/sigma on the diagonal.


    
%== GENERATE NOISY DATA ==================================================%
rng(0);
    % reset random number generator to make noise 
    % consistent between runs

if f_apx==1
    epsilon = sig.*randn(size(b0)); % noise vector (Gaussian approx.)
    b = sparse(b0+epsilon); % add noise to data
    b = max(round(b.*n_tot),0)./n_tot;
        % remove counts that would be below one and negative counts
else
    % Data vector generated using Poisson random numbers.
    % This will generally be similar to above.
    b = sig_gaus.*randn(size(b0))+...
        poissrnd(b0.*n_tot)./n_tot;
    b = max(b,0); % remove negative Gaussian noise
end


end

