
% MASS2RHO  Converts a mass-mobility distribution to an effective density-mobility distribution.
%  
%  Y = tools.mass2rho(X, GRID_X) converts a mass-mobility distribtuion to
%  an effective density-mobility distribution, using effective densities
%  that span from 100 -> 3,000 kg/m3 and 600 effective density elements.
%  
%  Y = tools.mass2rho(X, GRID_X, SPAN_RHO) exclitly states the range of
%  effective densities for converted space.
%  
%  Y = tools.mass2rho(..., N_RHO) explicitly states the number of elements
%  in the effective density dimension. 
%  
%  [Y, GRID_RHO] = tools.mass2rho(...) outputs the grid on which the
%  effective density-mobility distribution is defined.
%  
%  AUTHOR: Timothy Sipkens, 2019-05-17

function [y, grid_rho] = mass2rho(x, grid_x, span_rho, n_rho)

%-- Parse inputs -----------------------------------%
if ~exist('n_rho','var'); n_rho = []; end
if isempty(n_rho); n_rho = 600; end

if ~exist('span_rho','var'); span_rho = []; end
if isempty(span_rho); span_rho = [100,3000]; end
%---------------------------------------------------%


rho_min = span_rho(1);  % get span for effective density
rho_max = span_rho(2);
rho_n = logspace(log10(rho_min), ...
                 log10(rho_max), n_rho);  % discretize rho space

% Generate new grid for effective density-mobility.
grid_rho = Grid([rho_min,rho_max;grid_x.span(2,:)],...
    [n_rho,length(grid_x.edges{2})], 'logarithmic');

x_rs = grid_x.reshape(x);


%== Loop over mobility diameter ==========================================%
%   (i.e., consider conditional mass distributions)
n_d = grid_x.ne(2);
y = zeros(grid_x.ne(2), length(rho_n));

for ii=1:n_d
    
    % Initialize transformation kernel of appropriate size.
    c4 = zeros(grid_rho.ne(1),grid_x.ne(1));
    
    % Convert x nodes to effective density for iith mobility.
    rho_old = log10( ...
        6 .* grid_x.nodes{1} ./ (pi .* grid_x.edges{2}(ii) .^ 3) ...  % eq. for effective density
        .* 1e9);
    
    for jj=1:grid_x.ne(1)  % loop over masses
        c4(:,jj) = max(...
            min(log10(grid_rho.nodes{1}(2:end)), rho_old(jj+1)) - ... % upper bound
            max(log10(grid_rho.nodes{1}(1:(end-1))), rho_old(jj))... % lower bound
            ,0) ./ ...
            (log10(grid_rho.nodes{1}(2:end)) - ...
            log10(grid_rho.nodes{1}(1:(end-1)))); % normalize by rho bin size
    end
    
    % Multiply transformation by mass-mobility distr.
    y(ii,:) = c4 * x_rs(:,ii);
end
%=========================================================================%

% Convert to partial grid (due to missing entries).
% rho0 = 6/pi .* grid_x.span(1,1) .* 100.^(-3) .* 1e9;  % lower bound
% rho1 = 6/pi .* grid_x.span(1,2) .* 100.^(-3) .* 1e9;  % upper bound
% grid_rho = grid_rho.partial(...
%     log10([rho0, 100]), 3, ...
%     log10([rho1, 100]), 3);

% Format data for output
y = y';
y = y(:);
% y = grid_rho.full2partial(y);

end

