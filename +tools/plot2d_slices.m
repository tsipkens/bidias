
% PLOT2D_SLICES  Sweep through data using patch and creating a 3D plot.
% Author: Timothy Sipkens, 2020-04-06
%=========================================================================%

function [] = plot2d_slices(grid,x,cm,dim)

if ~exist('dim','var'); dim = []; end
if isempty(dim); dim = 1; end
    % dimension to sweep through
    % e.g. sweep through mass setpoints on standard grid, dim = 1

dim2 = setdiff([1,2],dim); % other dimension, dimension to plot


x_rs = grid.reshape(x); % reshape data
if dim==1; x_rs = x_rs'; end
min_x = max(log10(x))-3;

n1 = floor(size(cm,1)/grid.ne(dim));
n2 = length(cm)-grid.ne(dim)*n1+1;
cm2 = cm(n2:n1:end,:); % adjust colormap to appropriate size

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

view([-20,45,70]); % adjust view so slices are visible

end
