
% TWOMARK  Performs inversion using the iterative Twomey approach with intermediate smoothing.
% Author:  Timothy Sipkens, 2018-12-20
% Note:    Includes sub-functions responsible for smoothing.
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   Lb          Cholesky factorization of inverse covariance matrix
%   n           Length of first dimension of solution, used in smoothing
%   xi          Initial guess
%   iter        Max. number of iterations of Twomey_Markowski algorithm
%   opt_smooth  Type of smoothing to apply      (Optional, default is 'Buckley')
%   Sf          Smoothing parameter             (Optional, default is 1/300)
%   SIGMA_end   Mean square error for exit      (Optional, default is 1)
%
% Outputs:
%   x           Estimate
%=========================================================================%

function x = twomark(A,b,Lb,n,xi,iter,opt_smooth,Sf)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('Sf','var')
    Sf = 1/300; % matching Buckley et al.
elseif isempty(Sf)
    Sf = 1/300;
end

if ~exist('opt_smooth','var')
    opt_smooth = 'Buckley';
elseif isempty(opt_smooth)
    opt_smooth = 'Buckley';
end
%-------------------------------------------------------------------------%


iter_two = 150; % max number of iterations in Twomey pass
iter_2m = iter; % max number of iterations of Twomey_Markowski algorithm

x = xi;
x = invert.twomey(A,b,x,iter_two); % initial Towmey procedure
SIGMA = calc_mean_sq_error(Lb*A,x,Lb*b); % average square error for cases where b~= 0
R = roughness(x,n); % roughness vector

iter_two = 150; % max number of iterations in Twomey pass
for kk=1:iter_2m % iterate Twomey and smoothing procedure
    x_temp = x; % store temporarily for the case that roughness increases
    x = smooth_markowski(A,b,Lb,x,n,10,opt_smooth,Sf,SIGMA); % perform smoothing

    SIGMA_fun = @(x) calc_mean_sq_error(Lb*A,x,Lb*b);
    x = invert.twomey(A,b,x,iter_two,SIGMA_fun,SIGMA); % perform Twomey

    %-- Check roughness of solution ------------------%
    R(kk+1) = roughness(x,n);
    if R(kk+1)>1.03*R(kk) % exit if roughness has stopped decreasing
        disp(['Exited Twomey-Markowski loop after ',num2str(kk),...
            ' iterations because roughness increased by more than 3%.']);
        x = x_temp; % restore previous iteration
        break;
    end

    disp(['Completed iteration ',num2str(kk),' of the Twomey-Markowski loop.']);
    disp(' ');
end
disp('Completed Twomey-Markowski procedure.');
disp(' ');
disp(' ');

end
%=========================================================================%


%== SMOOTH_MARKOWSKI =====================================================%
%   Sub-funtion to perform smoothing required for Twomey-Markowski algorithm.
function [x,G_smooth,SIGMA] = smooth_markowski(A,b,Lb,x,n,iter,opt_smooth,Sf,SIGMA_end)

x_length = length(x);

if strcmp('Buckley',opt_smooth)
    G_smooth = G_Buckley(n,x_length,Sf);
else
    G_smooth = G_Markowski(n,x_length);
end

%-- Perform smoothing over multiple iterations ---------------------------%
jj = 1;
while jj<=iter
    x = G_smooth*x; % apply smoothing

    SIGMA = calc_mean_sq_error(Lb*A,x,Lb*b); % calculate mean square error
    if SIGMA>SIGMA_end % exit smoothing if mean square error exceeds end position
        disp(['SMOOTHING: Completed smoothing algorithm after ',num2str(jj),' iteration(s).']);
        return
    end

    jj = jj+1;
end

disp(['SMOOTHING algorithm did not necessarily converge after ',num2str(iter),' iteration(s).']);

end
%=========================================================================%


%== G_MARKOWSKI ==========================================================%
% Generates a smoothing matrix of the form based on original work of Markowski
function G_smooth = G_Markowski(n,x_length)

G_smooth = 0.5.*eye(x_length);
for jj=1:x_length
    if ~(mod(jj,n)==0)
        G_smooth(jj,jj+1) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end

    if ~(mod(jj-1,n)==0)
        G_smooth(jj,jj-1) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end

    if jj<=(x_length-n)
        G_smooth(jj,jj+n) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end

    if jj>n
        G_smooth(jj,jj-n) = 0.125;
    else
        G_smooth(jj,jj) = G_smooth(jj,jj)+0.125;
    end
end

end
%=========================================================================%


%== G_BUCKLEY ============================================================%
% Generates the smoothing matrix suggested by Buckley et al. (2017)
function G_smooth = G_Buckley(n,x_length,Sf)

norm_factor = 0.5+8*Sf;
G_smooth = 0.5/norm_factor.*eye(x_length);

for jj=1:x_length
    if mod(jj,n)==0 % right side
        if jj>(x_length-n) % bottom, right corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj-1) = 0.125;
            G_smooth(jj,jj-n) = 0.125;
        elseif jj<=n % top, right corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj-1) = 0.125;
            G_smooth(jj,jj+n) = 0.125;
        else % in the middle, right
            G_smooth(jj,jj-1) = 0.25;
            G_smooth(jj,jj) = 0.75;
        end
    elseif mod(jj-1,n)==0  % left side
        if jj>(x_length-n) % bottom, left corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj+1) = 0.125;
            G_smooth(jj,jj-n) = 0.125;
        elseif jj<=n % top, left corner
            G_smooth(jj,jj) = 0.75;
            G_smooth(jj,jj+1) = 0.125;
            G_smooth(jj,jj+n) = 0.125;
        else % in the middle, left
            G_smooth(jj,jj+1) = 0.25;
            G_smooth(jj,jj) = 0.75;
        end
    else % middle
        if jj>(x_length-n) % bottom, middle
            G_smooth(jj,jj-n) = 0.25;
            G_smooth(jj,jj) = 0.75;
        elseif jj<=n % top, middle
            G_smooth(jj,jj+n) = 0.25;
            G_smooth(jj,jj) = 0.75;
        else % in the middle, middle
            G_smooth(jj,jj-1) = Sf/norm_factor;
            G_smooth(jj,jj+1) = Sf/norm_factor;
            G_smooth(jj,jj+n) = Sf/norm_factor;
            G_smooth(jj,jj-n) = Sf/norm_factor;
            G_smooth(jj,jj+1+n) = Sf/norm_factor;
            G_smooth(jj,jj+1-n) = Sf/norm_factor;
            G_smooth(jj,jj-1+n) = Sf/norm_factor;
            G_smooth(jj,jj-1-n) = Sf/norm_factor;
        end
    end
end

end
%=========================================================================%


%== ROUGHNESS ============================================================%
%   Computes an estimate of the roughness of the solution.
%   This function is used for convergence in the Twomey-Markowski loop and
%   is based on the average, absolute value of the second derivative.
function R = roughness(x,n)

x = reshape(x,[n,length(x)/n]);
R = abs((x(3:end,:)+x(1:(end-2),:)-2.*x(2:(end-1),:)));
R = mean(mean(R));

end


%== CALC_MEAN_SQ_ERR =====================================================%
%   Calculates the mean squared error of the non-zero entries of b
function SIGMA = calc_mean_sq_error(A,x,b)

sqErr = (A*x-b).^2; % squared errors, normalized by expected error in b
SIGMA = mean(sqErr(b~=0)); % average square error for cases where b~= 0

end
%=========================================================================%
