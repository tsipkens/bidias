
% PLOT2D_SLICES  Sweep through data using patch and creating a 3D plot.
% Author: Timothy Sipkens, 2020-04-06
%=========================================================================%

function [] = plot2d_slices(grid,x,cm,dim)

% dimension to sweep through
% e.g. sweep through mass setpoints on standard grid, dim = 1
if ~exist('dim','var'); dim = []; end
if isempty(dim); dim = 1; end

% by default, use gray colormap if none specified
if ~exist('cm','var'); cm = []; end
if isempty(cm); cm = colormap('gray'); end

dim2 = setdiff([1,2],dim); % other dimension, dimension to plot


x_rs = grid.reshape(x); % reshape data
if dim==1; x_rs = x_rs'; end
min_x = max(log10(x))-3;

addpath('cmap'); % load cmap package to use `sweep_cmap(...)`
if isfile('cmap/cmap_sweep.m'); cm2 = cmap_sweep(grid.ne(dim), cm); % set color order to sweep through colormap
else; warning('The `cmap` package missing.'); end % if package is missing

clf;
plot3(log10(grid.edges{dim2}),... % plot data slices as patches
    log10(grid.edges{dim}(1)).*ones(1,grid.ne(dim2)),...
    max(log10(x_rs(:,1)'),min_x),'Color',cm2(1,:));
hold on;
for ii=2:grid.ne(dim)
    plot3(log10(grid.edges{dim2}),...
        log10(grid.edges{dim}(ii)).*ones(1,grid.ne(dim2)),...
        max(log10(x_rs(:,ii)'),min_x),'Color',cm2(ii,:));
end
hold off;
zlim([min_x,inf]);

view([-145,60]); % adjust view so slices are visible

end
