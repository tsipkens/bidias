function [y,rho_n,dNdlogrho] = ...
    massmob2rhomob(x,grid_x,n_rho)
%MASSMOB2RHOMOB Converts a mass-mobility distribution to a effective density-mobility distribution.
%   ...

%-- Parse inputs -----------------------------------%
if ~exist('n_rho','var') % if number of points in rho_n not specified
    n_rho = 100;
elseif isempty(n_rho)
    n_rho = 100;
end


rho_n = logspace(log10(100),log10(10000),n_rho);

x_rs = reshape(x,grid_x.ne);

drho = log(rho_n(2))-log(rho_n(1));
rho_n_low = exp(log(rho_n)-drho./2);
rho_n_high = exp(log(rho_n)+drho./2);

n_d = grid_x.ne(2);
y = zeros(grid_x.ne(2),length(rho_n));
for ii=1:n_d % loop over mobility diameter
    c3 = zeros(length(rho_n_high),grid_x.ne(1));
    for jj=1:length(rho_n_high)
        rho_old = 6.*grid_x.nodes{1}./(pi.*grid_x.edges{2}(ii).^3).*1e9;
        
        c0 = max(min(rho_old(2:end),rho_n_high(jj)),rho_n_low(jj));
        c1 = max(min(rho_old(1:(end-1)),rho_n_high(jj)),rho_n_low(jj));
        
        c2 = (c0-c1)./(rho_old(2:end)-rho_old(1:(end-1)));
        c3(jj,:) = c3(jj,:)+c2;
    end
    
    y(ii,:) = c3*x_rs(:,ii);
end

y = y';

dNdlogrho = sum(y,2);

end

