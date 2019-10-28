
[~,ind_max] = max(b_plot_rs,[],2);
d_max = grid_b.edges{2}(ind_max);
m_vec = grid_b.edges{1};

fun = @(x) (log10(m_vec) - log10(x(1)) - ...
    x(2).*log10(d_max./100)).^2;
x1 = lsqnonlin(fun,[1,2.2]);

hold on;
plot(log10(d_max),log10(grid_b.edges{1}),'ro');
plot(log10(grid_x.edges{2}),...
    log10(x1(1))+x1(2).*log10(grid_x.edges{2}./100),'r');
hold off;

m_100 = x1(1);
Dm = x1(2);
rho_100 = m_100/(pi*(100^3)/6)*1e9;
