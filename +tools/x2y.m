
% X2Y  A generic function to transform a distribution to a different space.
%  Examples include converting mass-mobility distributions to
%  effective-density mobility distributions or converting mrBC-mp
%  distributions to frBC-mp distributions. 
%  
%  Y = tools.x2y(X, GRID_X) converts a mass-mobility distribtuion to
%  an effective density-mobility distribution, using effective densities
%  that span from 100 -> 3,000 kg/m3 and 600 effective density elements.
%  
%  Y = tools.x2y(X, GRID_X, FUN, SPAN_Y) converts 
%  
%  Y = tools.x2y(..., N_Y) explicitly states the number of elements
%  in the effective density dimension. 
%  
%  [Y, GRID_Y] = tools.x2y(...) outputs the grid on which the
%  effective density-mobility distribution is defined.
%  
%  ------------------------------------------------------------------------
%  
%  NOTE: This function requires that FUN be monotonic within spans.
%  
%  AUTHOR: Timothy Sipkens, 2019-05-17

function [y, grid_y] = x2y(x, grid_x, fun, span_y, dim, n_y)

%-- Parse inputs -----------------------------------%
if ~exist('n_y', 'var'); n_y = []; end
if isempty(n_y); n_y = 600; end  % can be large as conversion is simple

if ~exist('span_y', 'var'); span_y = []; end
if isempty(span_y); span_y = [100, 3000]; end

if ~exist('fun', 'var'); fun = []; end
if isempty(fun)  % by default use mass-mobility distribution
    fun = @(a, b) 6 .* a ./ (pi .* b .^ 3) .* 1e9;
end

% By default, set first dimension as common dimension.
if ~exist('dim', 'var'); dim = []; end
if isempty(dim); dim = 2; end
dim2 = 3 - dim;  % other dimension, i.e., 2 -> 1 or 1 -> 2
%---------------------------------------------------%


y_min = span_y(1);  % get span for effective density
y_max = span_y(2);
y_n = logspace(log10(y_min), ...
               log10(y_max), n_y);  % discretize rho space

% Generate new grid for transformed space.
% Generate one new dimension while keeping dimension DIM.
if dim == 2
    grid_y = Grid([y_min,y_max; grid_x.span(2,:)],...
        [n_y, length(grid_x.edges{2})], 'logarithmic');
else
    grid_y = Grid([grid_x.span(1,:); y_min,y_max],...
        [length(grid_x.edges{1}), n_y], 'logarithmic');
end

x_rs = grid_x.reshape(x);


%== Loop over mobility diameter ==========================================%
%   (i.e., consider conditional mass distributions)
n_dim = grid_x.ne(dim);
y = zeros(grid_x.ne(dim), length(y_n));

for ii=1:n_dim
    
    % Initialize transformation kernel of appropriate size.
    T = zeros(grid_y.ne(dim2), grid_x.ne(dim2));
    
    % Convert x nodes to effective density for iith mobility.
    if dim == 2
        y_old = log10( ...
            fun(grid_x.nodes{1}, grid_x.edges{2}(ii)));
    else
        y_old = log10( ...
            fun(grid_x.edges{1}(ii), grid_x.nodes{2}));
    end
    
    % Reverse order if decreasing.
    % Assumes function in monotonic.
    if y_old(2) < y_old(1)
        y_old = fliplr(y_old);
        f_reverse = 1;
    else
        f_reverse = 0;
    end
    
    % Loop over masses computing overlap between new and old elements.
    for jj=1:grid_x.ne(dim2)
        T(:,jj) = max(...
            min(log10(grid_y.nodes{dim2}(2:end)), y_old(jj+1)) - ... % upper bound
            max(log10(grid_y.nodes{dim2}(1:(end-1))), y_old(jj))... % lower bound
            ,0) ./ ...
            (log10(grid_y.nodes{dim2}(2:end)) - ...
            log10(grid_y.nodes{dim2}(1:(end-1)))); % normalize by rho bin size
    end
    
    % Multiply transformation by mass-mobility distr.
    if f_reverse; T = fliplr(T); end
    if dim == 2
        y(ii,:) = T * x_rs(:,ii);
    else
        y(ii,:) = T * x_rs(ii,:)';
    end
end
%=========================================================================%


% Format data for output
if dim == 2; y = y'; end
y = y(:);

end

