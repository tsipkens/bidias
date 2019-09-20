
pha = Phantom.fit(grid_x,x_plot,'logn');

hold on;
plot(log10(grid_x.edges{2}),...
    log10(pha.p.m_100)+pha.p.Dm.*log10(grid_x.edges{2}./100),'w');
plot(log10(grid_x.edges{2}),...
    log10(pha.p.k)+...
    pha.p.Dm.*log10(grid_x.edges{2}),'g:');
hold off;
