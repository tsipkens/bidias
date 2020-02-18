
% FIT_MM_MODE Fit the mass-mobility relation to the mode of the SMPS scans.
% Author: Timothy Sipkens, 2019-01-21
%===========================================================%

function [p,d_max,m_vec] = fit_mm_mode(b,grid_b,grid_x,bool_plot)

if ~exist('bool_plot','var'); bool_plot = []; end
if isempty(bool_plot); bool_plot = 1; end


b_plot_rs = reshape(b,grid_b.ne);

%-- Preliminary analysis to find mode of SMPS scans --------%
[~,ind_max] = max(b_plot_rs,[],2);
m_vec = grid_b.edges{1};
d_max = grid_b.edges{2}(ind_max);


%-- Linear regression on median/mode diameter --------------%
opts = optimset('Display','none');
fun = @(x) log10(x(1)) + x(2).*log10(d_max./100);
min_fun = @(x) (fun(x)'-log10(m_vec)').^2;
x0 = [0.5,2.5];
x1 = lsqnonlin(min_fun,x0,[],[],opts);


%-- Plot results -------------------------------------------%
if bool_plot
    hold on;
    plot(log10(d_max),log10(m_vec),'c.');
    plot(log10(grid_x.edges{2}),...
        log10(x1(1)) + ...
        x1(2).*log10(grid_x.edges{2}./100),'c');
    hold off;
end


%-- Store morphological parameters -------------------------%
p.m_100 = x1(1);
p.Dm = x1(2);
p.rho_100 = p.m_100/(pi*(100^3)/6)*1e9;

end
