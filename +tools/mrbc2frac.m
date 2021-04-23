
% MRBC2FRAC  Converts a mp-mrBC distribution to a mp-frBC distribution.
%  
%  Y = tools.mrbc2frac(X, GRID_X) converts a mass-mobility distribtuion to
%  an effective density-mobility distribution, using effective densities
%  that span from 100 -> 3,000 kg/m3 and 600 effective density elements.
%  
%  Y = tools.mrbc2frac(X, GRID_X, SPAN_F) exclitly states the range of
%  mass fraction rBC for converted space.
%  
%  Y = tools.mrbc2frac(..., N_F) explicitly states the number of elements
%  in the mass fraction rBC dimension. 
%  
%  [Y, GRID_F] = tools.mrbc2frac(...) outputs the grid on which the
%  mp-frBC distribution is defined.
%  
%  AUTHOR: Timothy Sipkens, 2019-05-17

function [y, grid_f] = mrbc2frac(x, grid_x, span_f, n_f)

%-- Parse inputs -----------------------------------%
if ~exist('n_f','var'); n_f = []; end
if isempty(n_f); n_f = 600; end

if ~exist('span_f','var'); span_f = []; end
if isempty(span_f); span_f = [1e-3, 1]; end
%---------------------------------------------------%


f_min = span_f(1);  % get span for fraction rBC
f_max = span_f(2);
rho_n = logspace(log10(f_min), ...
                 log10(f_max), n_f); % discretize rho space

% Generate grid for mrBC-frBC.
grid_f = Grid([f_min, f_max; grid_x.span(2,:)],...
    [n_f, length(grid_x.edges{2})], 'logarithmic');

x_rs = grid_x.reshape(x);


%== Loop over mp =========================================================%
%   (i.e., consider conditional mrBC distributions)
n_mp = grid_x.ne(2);
y = zeros(grid_x.ne(2),length(rho_n));

for ii=1:n_mp  % loop over mobility diameter
    
    % Initialize transformation kernel of appropriate size.
    c4 = zeros(grid_f.ne(1),grid_x.ne(1));
    
    % convert x nodes to frBC for iith mp.
    f_old = log10(grid_x.nodes{1}./grid_x.edges{2}(ii));
    
    for jj=1:grid_x.ne(1)
        c4(:,jj) = max(...
            min(log10(grid_f.nodes{1}(2:end)),f_old(jj+1))-... % lower bound
            max(log10(grid_f.nodes{1}(1:(end-1))),f_old(jj))... % upper bound
            ,0)./...
            (log10(grid_f.nodes{1}(2:end))-log10(grid_f.nodes{1}(1:(end-1)))); % normalize by bin size
    end
    
    % Multiply transformation by mp-mrBC distr.
    y(ii,:) = c4 * x_rs(:,ii);
end
%=========================================================================%


% Format data for output
y = y';
y = y(:);

end

