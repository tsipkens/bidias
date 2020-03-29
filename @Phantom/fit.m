
% FIT  Fits a phantom to a given set of data, x, defined on a given grid. 
%      Outputs a fit phantom object.
% 
% Inputs:
%   x           Input data (2D distribution data)
%   vec_grid    A `Grid` object on which input data is defined or a vector
%               of the coordiantes of the elements
%   logr0       Initial guess for center of each mode, expressed as a vector:
%               logr0 = [dim1_mode1,dim2_mode1,dim1_mode2,dim2_mode2,...]
%
% Outputs: 
%   phantom     Output phantom object representing a fit bivariate lognormal
%   N           Scaling parameter that acts to scale the phantom to the data
%   y_out       Direct output from the fitting procedure
%   J           Jacobian of fitting procedure
%=============================================================%

function [phantom,N,y_out,J] = fit(x,vec_grid,logr0)

disp('[ Fitting phantom object... -------------]');

%-- Parse inputs ---------------------------------------------%
if isa(vec_grid,'Grid')
    grid = vec_grid;
    [~,vec1,vec2] = grid.vectorize();
else
    vec1 = vec_grid(:,1);
    vec2 = vec_grid(:,2);
    grid = [min(vec1),max(vec1);min(vec2),max(vec2)];
        % specify a span for a grid
end
if ~exist('logr0','var'); logr0 = []; end
%-------------------------------------------------------------%


corr2cov = @(sigma,R) diag(sigma)*R*diag(sigma);

fun_pha = @(y) y(1).*mvnpdf(log10([vec1,vec2]),[y(2),y(3)],...
    corr2cov([y(4),y(5)],[1,sin(y(6));sin(y(6)),1]));
y0 = [max(x),0,2.3,0.6,0.25,0.97];
    % [C,log10(mg),log10(dg),log10(sm),log10(sg),corr]

if ~isempty(logr0) % update centers of distributions, if specified
    y0(2:3) = logr0;
end

[y1,~,~,~,~,~,J] = lsqnonlin(@(y) fun_pha(y)-x, y0, ...
    [0,-10,-10,0,0,-1],[inf,10,10,10,3,1]);

mu = [y1(2),y1(3)];
sigma = [y1(4),y1(5)];
Sigma = corr2cov(sigma,[1,sin(y1(6));sin(y1(6)),1]);

phantom = Phantom('standard',grid,mu,Sigma);
phantom.type = 'standard-fit';

N = y1(1); % scaling parameter denoting total number of particles
y_out = y1;

disp(' ');
disp('[ Complete ------------------------------]');
disp(' ');
disp(' ');

end
