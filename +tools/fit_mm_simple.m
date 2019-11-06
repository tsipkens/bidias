
%-- Preliminary analysis to find mode of SMPS scans --------%
[~,ind_max] = max(b_plot_rs,[],2);
d_max = grid_b.edges{2}(ind_max);
d_max(1) = d_max(2); % account for noisy signals at the ends of the vector
d_max(end) = d_max(end-1);
m_vec = grid_b.edges{1};


%-- Fit normal distribution to data in log(d) space --------%
x2 = [];
opts = optimset('Display','none');
for bb=length(grid_b.edges{1}):-1:1
    fun = @(x) x(1).*normpdf(log10(grid_b.edges{2}),x(2),x(3));
    x0 = [max(b_plot_rs(bb,:))/5,log10(d_max(bb)),0.1];
    x2(bb,:) = lsqnonlin(@(x) (fun(x)-b_plot_rs(bb,:)).^2,x0,[],[],opts);
end
d_med = 10.^x2(:,2);

%-- Linear regression on median/mode diameter --------------%
fun = @(x) log10(x(1)) + x(2).*log10(d_med./100);
min_fun = @(x) (fun(x)-log10(m_vec)').^2;
x0 = [0.5,2.5];
x1 = lsqnonlin(min_fun,x0);%,[],[],opts);


%-- Plot results -------------------------------------------%
% figure(1); clf;
hold on;
plot(log10(d_med),log10(m_vec),'r.');
plot(log10(grid_x.edges{2}),...
    log10(x1(1)) + x1(2).*log10(grid_x.edges{2}./100),'r');
hold off;

m_100 = x1(1);
Dm = x1(2);
rho_100 = m_100/(pi*(100^3)/6)*1e9;
