
% PLOT2D_MARG   Plots x on the grid, with marginalized distributions.
% Author:       Timothy Sipkens, 2018-11-21
%=========================================================================%

function h = plot2d_marg(obj,x,obj_t,x_t)

subplot(4,4,[5,15]);
obj.plot2d(x);

x_m = obj.marginalize(x);


subplot(4,4,[1,3]);
marg_dim = 2;
stairs(log10(obj.nodes{marg_dim}),...
    [x_m{marg_dim},0],'k');
xlim(log10(obj.edges{2}([1,end])));

if nargin>2
    x_m_t = obj_t.marginalize(x_t);
    
    hold on;
    plot(log10(obj_t.nodes{marg_dim}),...
        [x_m_t{marg_dim},0],'color',[0.6,0.6,0.6]);
    hold off;
end


subplot(4,4,[8,16]);
marg_dim = 1;
stairs([0;x_m{marg_dim}],...
    log10(obj.nodes{marg_dim}),'k');
ylim(log10(obj.edges{1}([1,end])));

if nargin>2
    hold on;
    plot([0;x_m_t{marg_dim}],...
        log10(obj_t.nodes{marg_dim}),'color',[0.6,0.6,0.6]);
    hold off;
end


h = subplot(4,4,[5,15]);

end

