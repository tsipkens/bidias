
% FIT_MM_MEDIAN Fit the mass-mobility relation to the median.
% Author: Timothy Sipkens, 2019-01-21
%===========================================================%

function [p,d_med,m_vec] = fit_mm_median(b, grid_b, grid_x)

if ~exist('grid_x', 'var'); grid_x = []; end

% Infer whether to plot depending on if grid_x is supplied.
if isempty(grid_x); f_plot = 0;
else f_plot = 1; end


b_plot_rs = reshape(b,grid_b.ne);

%-- Preliminary analysis to find mode of SMPS scans --------%
%   Used as initial guess in lognormal distribution fitting
[~,ind_max] = max(b_plot_rs,[],2);
d_med = grid_b.edges{2}(ind_max);
d_med(1) = d_med(2); % account for noisy signals at the ends of the vector
d_med(end) = d_med(end-1);
m_vec = grid_b.edges{1};


%-- Fit normal distribution to data in log(d) space --------%
opts = optimset('Display','none');
for bb=length(grid_b.edges{1}):-1:1
    fun0 = @(x) x(1).*normpdf(log10(grid_b.edges{2}),x(2),x(3));
    x0 = [max(b_plot_rs(bb,:))/5,log10(d_med(bb)),0.1];
    [x2(bb,:),~,~,~,~,~,jacob] = ...
        lsqnonlin(@(x) (fun0(x)-b_plot_rs(bb,:)).^2,x0,[],[],opts);
    
    t0 = inv(jacob'*jacob);
    sx2(bb,1) = sqrt(t0(2,2)); % error in median, used to weight other analysis
end
d_med = 10.^x2(:,2);


%-- Linear regression on median/mode diameter --------------%
fun = @(x) log10(x(1)) + x(2).*log10(d_med./100);
min_fun = @(x) ((fun(x)-log10(m_vec)')./sx2).^2;
    % amounts to weighted least-squares
    
x0 = [0.5,2.5];
x1 = lsqnonlin(min_fun,x0,[],[],opts);


%-- Plot results -------------------------------------------%
if f_plot
    hold on;
    plot(log10(d_med),log10(m_vec),'r.');
    plot(log10(grid_x.edges{2}),...
        log10(x1(1)) + ...
        x1(2).*log10(grid_x.edges{2}./100),'r');
    hold off;
end


%-- Store morphological parameters -------------------------%
p.m_100 = x1(1);
p.Dm = x1(2);
p.rho_100 = p.m_100/(pi*(100^3)/6)*1e9;

end
