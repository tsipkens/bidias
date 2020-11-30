
% PLOT2D_PATCH  Sweep through data using patch and creating a 3D plot.
% Author: Timothy Sipkens, 2020-04-02
%=========================================================================%

function [] = plot2d_patch(grid, x, cm, dim, opts)

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
cm2 = cm(n2:n1:end,:);  % adjust colormap to appropriate size


plotter = @patch;  % use `patch` by default

% If plotting lines instead of patches.
if opts.f_line==1
    plotter = @(x,y,z,cm) plot3(x,y,z,'Color',cm);
end


% Clear figure and 
% proceed to loop through slices.
clf;
hold on;
for ii=1:grid.ne(dim)
    p = plotter(grid.edges{dim2}([1,1:end,end]), ...
        grid.edges{dim}(ii) .* ones(1,grid.ne(dim2)+2), ...
        [min_x, max(x_rs(:,ii)', min_x), min_x], ...
        cm2(ii,:));
end
hold off;
zlim([min_x,inf]);

set(gca, 'ZScale', 'log', 'YScale', 'log', 'XScale', 'log');

view([-145,60]); % adjust view so slices are visible

end
