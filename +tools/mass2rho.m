
% MASS2RHO  Converts a mass-mobility distribution to a effective density-mobility distribution.
% Author:   Timothy Sipkens, 2019-05-17
%=========================================================================%

function [y,grid_rho] = mass2rho(x,grid_x,span_rho,n_rho)

%-- Parse inputs -----------------------------------%
if ~exist('n_rho','var'); n_rho = []; end
if isempty(n_rho); n_rho = 600; end

if ~exist('span_rho','var'); span_rho = []; end
if isempty(span_rho); span_rho = [10^1.5,10000]; end
%---------------------------------------------------%

rho_min = span_rho(1);
rho_max = span_rho(2);
rho_n = logspace(log10(rho_min),log10(rho_max),n_rho);
grid_rho = Grid([rho_min,rho_max;grid_x.span(2,:)],...
    [n_rho,length(grid_x.edges{2})],'logarithmic'); % should be uniform basis

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
y = y(:);

end

