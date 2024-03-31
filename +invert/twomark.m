
% TWOMARK  Performs inversion using the iterative Twomey-Markowski approach.
%  This add an intermediate smoothing step to the standard Twomey routine. 
%  
%  X = invert.twomark(A, B, LB, N, XI, ITER) uses the model matrix, A;
%  data, B; Cholesky factorization of inverse covariance matrix, LB; 
%  length of first dimension of solution, used in smoothing, N;  
%  initial guess, XI; and maximum number of iterations of the algorithm,
%  ITER, udring inversion. 
%  
%  X = invert.twomark(A, B, LB, GRID, XI, ITER) replaces N with a Grid used
%  to inform on smoothing. 
%  
%  X =invert.twomark(..., OPT_SMOOTH, SF) adds optional parameters:
%  OPT_SMOOTH, the type of smoothing (default value: 'Buckley'), and 
%  SF, a smoothing parameter (default value: 1/300).
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR:  Timothy Sipkens, 2018-12-20

function x = twomark(A, b, Lb, n_grid, xi, iter, opt_smooth, Sf)

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
    if isa(n_grid, 'PartialGrid')
        opt_smooth = 'Grid';
    elseif ~strcmp(opt_smooth, 'Grid')
        n_grid = n_grid.ne(1);
    end
end
%-------------------------------------------------------------------------%


iter_two = 150; % max number of iterations in Twomey pass
iter_mh = iter; % max number of iterations of Twomey_Markowski algorithm

x = xi;
x = invert.twomey(A,b,x,iter_two); % initial Towmey procedure
SIGMA = calc_mean_sq_error(Lb*A,x,Lb*b); % average square error for cases where b~= 0
R = roughness(x, n_grid); % roughness vector

iter_two = 150; % max number of iterations in Twomey pass
for kk=1:iter_mh % iterate Twomey and smoothing procedure
    x_temp = x; % store temporarily for the case that roughness increases
    x = invert.markowski(A,b,Lb,x,n_grid,10,opt_smooth,Sf,SIGMA); % perform smoothing

    SIGMA_fun = @(x) calc_mean_sq_error(Lb*A,x,Lb*b);
    x = invert.twomey(A,b,x,iter_two,SIGMA_fun,SIGMA); % perform Twomey

    %-- Check roughness of solution ------------------%
    R(kk+1) = roughness(x,n_grid);
    if R(kk+1)>1.03*R(kk) % exit if roughness has stopped decreasing
        disp(['Exited Twomey-Markowski loop after ',num2str(kk),...
            ' iterations because roughness increased by more than 3%.']);
        x = x_temp; % restore previous iteration
        break;
    end
end
disp('Completed Twomey-Markowski procedure:');
disp(['iter = ',num2str(iter_mh)]);
disp(' ');

end
%=========================================================================%



%== ROUGHNESS ============================================================%
%   Computes an estimate of the roughness of the solution.
%   This function is used for convergence in the Twomey-Markowski loop and
%   is based on the average, absolute value of the second derivative.
function R = roughness(x,n_grid)

if isa(n_grid,'Grid')
    x = n_grid.reshape(x);
else
    x = reshape(x,[n_grid, length(x)/n_grid]);
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
