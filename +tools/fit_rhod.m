
% FIT_RHOD Fit the mass-mobility relation to the modes, plus plot other fits.
% Performs a mode-based analysis based on the effective density as well as 
% plotting mobility distribution-based fit, and the universal expression.
% Author: Timothy Sipkens, 2019-01-21
%=========================================================================%

function [p,d_max,m_vec] = fit_rhod(b,grid_b,grid_x,pha,bool_plot)

if ~exist('bool_plot','var'); bool_plot = []; end
if isempty(bool_plot); bool_plot = 1; end


b_plot_rs = reshape(b,grid_b.ne);

%-- Preliminary analysis to find mode of SMPS scans --------%
[~,ind_max] = max(b_plot_rs,[],2);
m_vec = grid_b.edges{1};
d_max = grid_b.edges{2}(ind_max);
rho_max = 1e9.*m_vec./(pi/6.*(d_max).^3);
    % mode-based analysis, with transformed data


fun = @(x) log10(rho_max) - log10(x(1)) - ...
    (x(2)-3).*log10(d_max./100);
x1 = lsqnonlin(fun,[500,0.1]);
    % fit to transformed data
    % should be similar to prev. mode-based analysis


if bool_plot
    hold on;
    plot(log10(d_max),log10(rho_max),'c.');
    plot(log10(grid_x.edges{2}),log10(x1(1))+...
        (x1(2)-3).*log10(grid_x.edges{2}./100),'c');
        % mode-based
    
    plot(log10(grid_x.edges{2}),...
        log10(1e9.*pha.p.m_100.*(grid_x.edges{2}./100).^pha.p.Dm./...
        (pi/6.*grid_x.edges{2}.^3)),'w');
        % distribution-based
    
    plot(log10(grid_x.edges{2}),...
        log10(510.*(grid_x.edges{2}./100).^-0.52),'g');
        % universal expression
    hold off;
end


p.rho_100 = x1(1);
p.Dm = x1(2);

end
