
% MRBC2FRAC  Converts a mrBC-mp distribution to a fraction mrBC-mp distribution.
% Author:   Timothy Sipkens, 2019-05-17
%=========================================================================%

function [y,grid_f] = mrbc2frac(x,grid_x,span_f,n_f)

%-- Parse inputs -----------------------------------%
if ~exist('n_f','var'); n_f = []; end
if isempty(n_f); n_f = 600; end

if ~exist('span_f','var'); span_f = []; end
if isempty(span_f); span_f = [1e-3,1]; end
%---------------------------------------------------%


f_min = span_f(1); % get span for fraction rBC
f_max = span_f(2);
rho_n = logspace(log10(f_min),log10(f_max),n_f); % discretize rho space
grid_f = Grid([f_min,f_max;grid_x.span(2,:)],...
    [n_f,length(grid_x.edges{2})],'logarithmic');
    % generate grid

x_rs = grid_x.reshape(x);

n_d = grid_x.ne(2);
y = zeros(grid_x.ne(2),length(rho_n));
for ii=1:n_d % loop over mobility diameter
    c4 = zeros(grid_f.ne(1),grid_x.ne(1));
    f_old = log10(grid_x.nodes{1}./grid_x.edges{2}(ii));
            % convert x nodes to effective density for iith mobility
    
    for jj=1:grid_x.ne(1)
        c4(:,jj) = max(...
            min(log10(grid_f.nodes{1}(2:end)),f_old(jj+1))-... % lower bound
            max(log10(grid_f.nodes{1}(1:(end-1))),f_old(jj))... % upper bound
            ,0)./...
            (log10(grid_f.nodes{1}(2:end))-log10(grid_f.nodes{1}(1:(end-1)))); % normalize by rho_old bin size
    end
    
    y(ii,:) = c4*x_rs(:,ii);
end

y = y';
y = y(:);

end

