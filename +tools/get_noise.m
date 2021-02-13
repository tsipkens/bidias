
% GET_NOISE Simulates Poisson-Gaussian noise in signals.
%  Modeled after Sipkens et al., Appl. Opt. (2017).
%  DOI: https://doi.org/10.1364/AO.56.008436
%  
%  B = tools.get_noise(B0) adds Poisson-Gaussian noise to the raw signal
%  B0 (representing raw counts) and outputs data corrupted wiht noise to B.
%  
%  B = tools.get_noise(B0,N_TOT) adds the option to use scaled data as
%  input. By default, N_TOT = 1, which means that B0 is counts. For
%  example, this usage can corrupt a PDF specified by B0 with
%  Poisson-Gaussian noise if the total counts over the PDF was N_TOT. 
%  Note that tools.get_noise(N_TOT.*B0)./N_TOT and 
%  tools.get_noise(B0,N_TOT) will give the same result. 
%  
%  B = tools.get_noise(B0,N_TOT,GAM0) adds an input to specify the level of
%  Gaussian background noise. The quantity is phrased relative to the
%  maximum value of N_TOT.*B0. By default, GAM0 = 1e-4, which corresponds 
%  to a Gaussian noise level of 0.01% of the the peak noise.
%  
%  B = tools.get_noise(...,F_APX) adds a flag for whether to use a Gaussian 
%  approximation is used for the Poisson noise. 
%  
%  [B,LB] = tools.get_noise(...) adds an output for the Cholesky 
%  factorization of the inverse of the covariance matrix, assuming the 
%  Poisson noise can be approximated as Gaussian. 
%  
%  ------------------------------------------------------------------------
%  
%  NOTE: Note that Lb is the Cholesky factorization of the
%  inverse of the covariance matrix for the Gaussian
%  approximation. Though, this covariance will 
%  well-approximates the covariance information for the
%  true Poisson case. 
%  
%  AUTHOR: Timothy Sipkens, 2019-12-08

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

if f_apx==1 % if using Gaussian approximation
    epsilon = sig.*randn(size(b0)); % noise vector (Gaussian approx.)
    b = sparse(b0+epsilon); % add noise to data
    b = max(round(b.*n_tot),0)./n_tot;
        % remove counts that would be below one and negative counts

else % if using Poisson distribution directly
    % Data vector generated using Poisson random numbers.
    % This will generally be similar to above.
    b = sig_gaus.*randn(size(b0))+...
        poissrnd(b0.*n_tot)./n_tot;
    b = max(b,0); % remove negative Gaussian noise
end


end

