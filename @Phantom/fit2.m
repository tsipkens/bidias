
% FIT2  Fits a multimodal phantom to a given set of data, x, defined on a given grid. 
%       Outputs a fit phantom object.
% 
% Inputs:
%   x         Input data (2D distribution data)
%   vec_grid  A `Grid` object on which input data is defined or a vector
%             of the coordiantes of the elements
%   n_modes   Number of modes to fit to the data.
%
% Outputs: 
%   phantom   Output phantom object representing a fit bivariate lognormal.
%   N         Scaling parameter that acts to scale the phantom to the data.
%=============================================================%

function [phantom,N] = fit2(x,vec_grid,n_modes,logr0)

disp('Fitting phantom object...');

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

y0 = [];
y1 = [];
yup = [];
ylow = [];
for ii=1:n_modes
    y0 = [y0,log(max(x))-3*(ii-1)/n_modes-3.5,...
        3+ii/n_modes,1.8+ii/n_modes,0.15,0.15,-0.9];
    ylow = [ylow,-inf,-3,-3,0,0,-pi];
    yup = [yup,inf,5,5,0.2,0.3,pi];
        % [C,log10(mg),log10(dg),log10(sm),log10(sg),corr]
end
if ~isempty(logr0) % update centers of distributions, if specified
    y0(2:6:end) = logr0(1:2:end);
    y0(3:6:end) = logr0(2:2:end);
end
for ii=1:n_modes % prior std. dev.
    y1 = [y1,inf,1e-2.*y0(6*(ii-1)+[2,3]),y0(6*(ii-1)+[4,5]),10];
        % [C,log10(mg),log10(dg),log10(sm),log10(sg),corr]
end

% opts = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1e4,'MaxIterations',1e3);
% y2 = lsqnonlin(@(y) fun_pha(y,vec1,vec2,n_modes,corr2cov)-...
%     x, y0, ylow, yup);
y2 = lsqnonlin(@(y) [max(log(fun_pha(y,vec1,vec2,n_modes,corr2cov)),max(log(x)-4))-...
    max(log(x),max(log(x)-4));...
    (y-y0)'./y1'], y0, ylow, yup);

for ii=0:(n_modes-1)
    mu{ii+1} = [y2(6*ii+2),y2(6*ii+3)];
    sigma{ii+1} = [y2(6*ii+4),y2(6*ii+5)];
    Sigma{ii+1} = corr2cov(sigma{ii+1},[1,sin(y2(6*ii+6));sin(y2(6*ii+6)),1]);
end

phantom = Phantom('standard',grid,mu,Sigma);
phantom.type = 'standard-fit';

N = exp(y2(1:6:end)); % scaling parameter denoting total number of particles
disp('Complete.');
disp(' ');

end


function x = fun_pha(y,vec1,vec2,n_modes,corr2cov)

x = 0;
for ii=0:(n_modes-1)
    x = x + exp(y(6*ii+1)).*mvnpdf(log10([vec1,vec2]),[y(6*ii+2),y(6*ii+3)],...
        corr2cov([y(6*ii+4),y(6*ii+5)],[1,sin(y(6*ii+6));sin(y(6*ii+6)),1]));
            % cos() correlation is to better span -1 to 1
end

end

