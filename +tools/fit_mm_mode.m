
% FIT_MM_MODE Fit the mass-mobility relation to the mode of the SMPS scans.
% Author: Timothy Sipkens, 2019-01-21
%===========================================================%

function [p,d_max,m_vec] = fit_mm_mode(b, grid_b, grid_x)

if ~exist('grid_x', 'var'); grid_x = []; end

% Infer whether to plot depending on if grid_x is supplied.
if isempty(grid_x); f_plot = 0;
else f_plot = 1; end


b_plot_rs = reshape(b, grid_b.ne);

%-- Preliminary analysis to find mode of SMPS scans --------%
[b_max, ind_max] = max(b_plot_rs, [], 2);

m_vec = grid_b.edges{1};  % form [m,d] pairs for peaks
d_max = grid_b.edges{2}(ind_max);

% Ignore peaks below 5% of max.
ind_rmv = b_max < (0.02 .* max(b_max));
d_max(ind_rmv) = [];
m_vec(ind_rmv) = [];



%-- Linear regression on median/mode diameter --------------%
fun = @(x) log10(x(1)) + x(2).*log10(d_max./100);  % function for line through results
min_fun = @(x) fun(x)' - log10(m_vec)';  % function to optimize

x0 = [0.5, 2.5];  % starting point, Dm = 2.5
opts = optimset('Display','none');
x1 = lsqnonlin(min_fun,x0,[],[],opts);  % fit a line


%-- Plot results -------------------------------------------%
if f_plot
    grid_b.plot2d_scatter(b);
    hold on;
    plot(d_max, m_vec, 'c.');
    plot(grid_x.edges{2}, ...
        10 .^ (log10(x1(1)) + ...
        x1(2).*log10(grid_x.edges{2}./100)), 'c');
    hold off;
end


%-- Store morphological parameters -------------------------%
p.m_100 = x1(1);
p.Dm = x1(2);
p.rho_100 = p.m_100/(pi*(100^3)/6)*1e9;

end
