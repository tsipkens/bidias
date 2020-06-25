
% MASS2RHO  Converts a mass-mobility distribution to an effective density-mobility distribution.
% Author:   Timothy Sipkens, 2019-05-17
%=========================================================================%

function [y,grid_rho] = mass2rho(x,grid_x,span_rho,n_rho)

%-- Parse inputs -----------------------------------%
if ~exist('n_rho','var'); n_rho = []; end
if isempty(n_rho); n_rho = 600; end

if ~exist('span_rho','var'); span_rho = []; end
if isempty(span_rho); span_rho = [100,3000]; end
%---------------------------------------------------%


rho_min = span_rho(1); % get span for effective density
rho_max = span_rho(2);
rho_n = logspace(log10(rho_min),log10(rho_max),n_rho); % discretize rho space
grid_rho = Grid([rho_min,rho_max;grid_x.span(2,:)],...
    [n_rho,length(grid_x.edges{2})],'logarithmic');
    % generate grid

x_rs = grid_x.reshape(x);

n_d = grid_x.ne(2);
y = zeros(grid_x.ne(2),length(rho_n));
for ii=1:n_d % loop over mobility diameter (consider conditional mass distributions)
    c4 = zeros(grid_rho.ne(1),grid_x.ne(1)); % initialize transformation kernel
    rho_old = log10(6.*grid_x.nodes{1}./(pi.*grid_x.edges{2}(ii).^3).*1e9);
            % convert x nodes to effective density for iith mobility
    
    for jj=1:grid_x.ne(1)
        c4(:,jj) = max(...
            min(log10(grid_rho.nodes{1}(2:end)),rho_old(jj+1))-... % upper bound
            max(log10(grid_rho.nodes{1}(1:(end-1))),rho_old(jj))... % lower bound
            ,0)./...
            (log10(grid_rho.nodes{1}(2:end))-log10(grid_rho.nodes{1}(1:(end-1)))); % normalize by rho bin size
    end
    
    y(ii,:) = c4*x_rs(:,ii);
end

y = y';
y = y(:);

end

