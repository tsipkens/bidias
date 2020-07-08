
% PLOT2D_PATCH  Sweep through data using patch and creating a 3D plot.
% Author: Timothy Sipkens, 2020-04-02
%=========================================================================%

function [] = plot2d_patch(grid,x,cm,dim)

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

n1 = floor(size(cm,1)/grid.ne(dim));
n2 = length(cm)-grid.ne(dim)*n1+1;
cm2 = cm(n2:n1:end,:); % adjust colormap to appropriate size

clf;
p = patch(log10(grid.edges{dim2}([1,1:end,end])),... % plot data slices as patches
    log10(grid.edges{dim}(1)).*ones(1,grid.ne(dim2)+2),...
    [min_x,max(log10(x_rs(:,1)'),min_x),min_x],cm2(1,:));
p.FaceAlpha = 1;
hold on;
for ii=2:grid.ne(dim)
    p = patch(log10(grid.edges{dim2}([1,1:end,end])),...
        log10(grid.edges{dim}(ii)).*ones(1,grid.ne(dim2)+2),...
        [min_x,max(log10(x_rs(:,ii)'),min_x),min_x],cm2(ii,:));
    p.FaceAlpha = 1;%1-ii/(grid.ne(dim));
end
hold off;
zlim([min_x,inf]);

view([-145,60]); % adjust view so slices are visible

end
