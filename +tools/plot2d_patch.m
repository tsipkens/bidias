
% PLOT2D_PATCH  Sweep through data using patch and creating a 3D plot.
% Author: Timothy Sipkens, 2020-04-02
%=========================================================================%

function [] = plot2d_patch(grid,x,cm)

x_rs = grid.reshape(x);
min_x = max(log10(x))-3;

n1 = ceil(grid.ne(1)./20);
n2 = floor(grid.ne(1)/n1);
n3 = floor(length(cm)/n2);
cm2 = cm(1:n3:end,:);

patch(log10(grid.edges{1}),...
    log10(grid.edges{2}(1)).*ones(1,grid.ne(1)),...
    max(log10(x_rs(:,1)'),min_x),cm2(1,:));
for ii=2:grid.ne(2)
    hold on;
    patch(log10(grid.edges{1}),...
        log10(grid.edges{2}(ii)).*ones(1,grid.ne(1)),...
        max(log10(x_rs(:,ii)'),min_x),cm2(ii,:));
    hold off;
end

view([-20,45,70]);

end
