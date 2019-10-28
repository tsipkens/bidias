
% [~,ind_max] = max(b_plot_rs,[],2);
d_max = grid_b.edges{2}(ind_max);
rho_max = 1e9.*m_vec./(pi/6.*(d_max).^3);

fun = @(x) log10(rho_max) - log10(x(1)) - ...
    (x(2)-3).*log10(d_max./100);
x1 = lsqnonlin(fun,[500,0.1]);

hold on;
plot(log10(d_max),log10(rho_max),'ro');
plot(log10(grid_x.edges{2}),log10(x1(1))+...
    (x1(2)-3).*log10(grid_x.edges{2}./100),'r');
plot(log10(grid_x.edges{2}),...
    log10(1e9.*pha.p.m_100.*(grid_x.edges{2}./100).^pha.p.Dm./...
    (pi/6.*grid_x.edges{2}.^3)),'w');
plot(log10(grid_x.edges{2}),...
    log10(510.*(grid_x.edges{2}./100).^-0.52),'y');
hold off;

rho_100_0 = x1(1);
Dm_0 = x1(2);
