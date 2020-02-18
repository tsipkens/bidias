
% PLOT2D_MARG_B Plots x as a 2D function on the grid, with marginalized distributions.
% Author: Timothy Sipkens, 2020-02-03
%==================================================================%

function [h,x_m] = plot2d_marg_b(grid,x,obj_t,x_t)

subplot(4,6,[8,23]);
grid.plot2d(x);

x_m = grid.marginalize(x);

%-- Plot marginal distribution (dim 2) -----------------------%
subplot(4,6,[2,5]);
marg_dim = 2;
stairs(grid.nodes{marg_dim},...
    [x_m{marg_dim},0],'k');
xlim([min(grid.edges{marg_dim}),max(grid.edges{marg_dim})]);
set(gca,'XScale','log');

if nargin>2 % also plot marginal of the true distribution
    x_m_t = obj_t.marginalize(x_t);

    hold on;
    plot(obj_t.nodes{marg_dim},...
        [x_m_t{marg_dim},0],'color',[0.6,0.6,0.6]);
    hold off;
end

%-- Plot marginal distribution (dim 1) -----------------------%
subplot(4,6,[12,24]);
marg_dim = 1;
stairs([0;x_m{marg_dim}],...
    grid.nodes{marg_dim},'k');
ylim([min(grid.edges{marg_dim}),max(grid.edges{marg_dim})]);
set(gca,'YScale','log');

if nargin>2 % also plot marginal of the true distribution
    hold on;
    plot([0;x_m_t{marg_dim}],...
        obj_t.nodes{marg_dim},'color',[0.6,0.6,0.6]);
    hold off;
end

subplot(4,6,[8,23]);
if nargout>0; h = gca; end

colorbar('position',[0.1,0.1,0.05,0.6]);

end


