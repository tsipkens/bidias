
% PLOT2D_STACK  Sweep through data and stack results.
% Author: Timothy Sipkens, 2020-04-02
%=========================================================================%

function [] = plot2d_stack(grid, x, cm, dim)

%-- Parse inputs ---------------------------------------------------------%
% dimension to sweep through
% e.g. sweep through mass setpoints on standard grid, dim = 1
if ~exist('dim','var'); dim = []; end
if isempty(dim); dim = 1; end

% by default, use gray colormap if none specified
if ~exist('cm','var'); cm = []; end
if isempty(cm); cm = colormap('gray'); end

dim2 = setdiff([1,2],dim); % other dimension, dimension to plot

if ~exist('opts', 'var'); opts = struct(); end
if ~isfield(opts, 'f_line'); opts.f_line = 0; end  % by default, apply patch
%-------------------------------------------------------------------------%


x_rs = grid.reshape(x);  % reshape data
if dim==1; x_rs = x_rs'; end
% min_x = 10.^(max(log10(x))-4);  % data minimum, applied due to log scale
min_x = 0;
dsc = mean(log(grid.edges{dim}(2:end)) - log(grid.edges{dim}(1:(end-1))));
sc = (exp(log(grid.edges{dim}) + dsc) - grid.edges{dim}) ./ max(x).* 0.9;


% Set up colormap for sweep.
n1 = floor(size(cm,1)/grid.ne(dim));
n2 = length(cm)-grid.ne(dim)*n1+1;
cm2 = cm(n2:n1:end,:); % adjust colormap to appropriate size

cla;
hold on;
for ii=1:grid.ne(dim)   
    p = tools.stairs(grid.edges{dim2}(1:end), ...
        sc(ii) .* x_rs(:,ii)' + grid.edges{dim}(ii), ...
        0, '-', 'Color', cm2(ii,:));
    plot(grid.edges{dim2}(1:end), ...
        (min_x + grid.edges{dim}(ii)) .* ones(grid.ne(dim2), 1), ...
        ':', 'Color', cm2(ii,:));
    text(grid.edges{dim2}(1) .* 1.05, ...
        min_x + grid.edges{dim}(ii), ...
        num2str(grid.edges{dim}(ii)), ...
        'Color', cm2(ii,:), ...
        'VerticalAlignment', 'bottom');
end
hold off;

set(gca, 'YScale', 'log', 'XScale', 'log');
xlim(grid.span(dim2,:));

end
