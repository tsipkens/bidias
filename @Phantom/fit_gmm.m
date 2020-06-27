
% FIT_GMM  Fits a phantom to a given set of data, x, defined on a given grid.
%   Uses sampling and k-means to fit a Guassian mixture model.
%   Outputs a fit phantom object.
% 
% Inputs:
%   x           Input data (2D distribution data)
%   grid        A `Grid` object on which input data is defined
%   k           Number of modes
%
% Outputs: 
%   phantom     Output phantom object representing a fit bivariate lognormal
%   N           Scaling parameter that acts to scale the phantom to the data
%   y_out       Direct output from the fitting procedure
%   J           Jacobian of fitting procedure
%=============================================================%

function [phantom,N,s] = fit_gmm(x,grid_vec,k)

disp(' ');
disp('[ Fitting phantom object... ===============]');
disp(' ');

disp('Sampling from distribution...');
if isa(grid_vec,'Grid') % get samples from distribution
    s = tools.sample_dist_a(x,grid_vec);
    [~,tot] = grid_vec.marginalize(x);
else
    s = tools.sample_dist_c(x,grid_vec);
    tot = sum(x);
end
disp('Sampling complete.');
disp(' ');

disp('Estimating a Gaussian mixture model...');
GMModel = fitgmdist(s,k);
N = GMModel.ComponentProportion.*tot;
disp('Complete');
disp(' ');

phantom = Phantom('standard',grid_vec,GMModel.mu,...
    GMModel.Sigma,GMModel.ComponentProportion);


disp('[ Complete ================================]');
disp(' ');

end

