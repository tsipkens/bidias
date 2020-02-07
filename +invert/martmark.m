
% MARTMARK  Performs inversion using the iterative MART approach with intermediate smoothing.
% Author: Arash Naseri, Timothy Sipkens, 2020-02-06
%-------------------------------------------------------------------------%
% Inputs:
%   A           Model matrix
%   b           Data
%   Lb          Cholesky factorization of inverse covariance matrix
%   n           Length of first dimension of solution, used in smoothing
%   xi          Initial guess
%   iter        Max. number of iterations of MART_Markowski algorithm
%   opt_smooth  Type of smoothing to apply      (Optional, default is 'Buckley')
%   Sf          Smoothing parameter             (Optional, default is 1/300)
%   SIGMA_end   Mean square error for exit      (Optional, default is 1)
%
% Outputs:
%   x           Estimate
%=========================================================================%

function x = martmark(A,b,Lb,n_grid,xi,iter,opt_smooth,Sf)


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

if isa(n_grid,'Grid') % handling of grid input
    if n_grid.ispartial==1
        opt_smooth = 'Grid';
    elseif ~strcmp(opt_smooth,'Grid')
        n_grid = n_grid.ne(1);
    end
end
%-------------------------------------------------------------------------%


iter_two = 150; % max number of iterations in MART pass
iter_2m = iter; % max number of iterations of MART-Markowski algorithm

x = xi;
x = invert.twomey(A,b,x,iter_two); % initial Towmey procedure
SIGMA = calc_mean_sq_error(Lb*A,x,Lb*b); % average square error for cases where b~= 0
R = roughness(x,n_grid); % roughness vector

iter_two = 150; % max number of iterations in MART pass
for kk=1:iter_2m % iterate MART and smoothing procedure
    x_temp = x; % store temporarily for the case that roughness increases
    x = invert.markowski(A,b,Lb,x,n_grid,10,opt_smooth,Sf,SIGMA); % perform smoothing

    SIGMA_fun = @(x) calc_mean_sq_error(Lb*A,x,Lb*b);
    x = invert.mart(A,b,x,iter_two,SIGMA_fun,SIGMA); % perform MART

    %-- Check roughness of solution ------------------%
    R(kk+1) = roughness(x,n_grid);
    if R(kk+1)>1.03*R(kk) % exit if roughness has stopped decreasing
        disp(['Exited MART-Markowski loop after ',num2str(kk),...
            ' iterations because roughness increased by more than 3%.']);
        x = x_temp; % restore previous iteration
        break;
    end

    disp(['Completed iteration ',num2str(kk),' of the MART-Markowski loop.']);
    disp(' ');
end
disp('Completed MART-Markowski procedure.');
disp(' ');
disp(' ');

end
%=========================================================================%



%== ROUGHNESS ============================================================%
%   Computes an estimate of the roughness of the solution.
%   This function is used for convergence in the MART-Markowski loop and
%   is based on the average, absolute value of the second derivative.
function R = roughness(x,n_grid)

if isa(n_grid,'Grid')
    x = n_grid.reshape(x);
else
    x = reshape(x,[n_grid,length(x)/n_grid]);
end

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
