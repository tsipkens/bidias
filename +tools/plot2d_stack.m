
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
min_x = 10.^(max(log10(x))-4);  % data minimum, applied due to log scale


% Set up colormap for sweep.
n1 = floor(size(cm,1)/grid.ne(dim));
n2 = length(cm)-grid.ne(dim)*n1+1;
cm2 = cm(n2:n1:end,:); % adjust colormap to appropriate size

clf;
hold on;
for ii=1:grid.ne(dim)
    p = patch(grid.edges{dim2}([1,1:end,end]), ...
        10.^(log10([min_x, max(x_rs(:,ii)', min_x), min_x]) + 8.*log10(grid.edges{dim}(ii))), ...
        cm2(ii,:), ...
        'EdgeColor', 'w');
    plot(grid.edges{dim2}(1:end), ...
        10.^(log10(min_x) + 8.*log10(grid.edges{dim}(ii))) .* ones(grid.ne(dim2), 1), ...
        'Color', cm2(ii,:));
    text(grid.edges{dim2}(1) .* 1.07, ...
        10.^(log10(min_x) + 8.*log10(grid.edges{dim}(ii))), ...
        num2str(grid.edges{dim}(ii)), ...
        'Color', cm2(ii,:), ...
        'VerticalAlignment', 'bottom');
end
hold off;

set(gca, 'YScale', 'log', 'XScale', 'log');

end
