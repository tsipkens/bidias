
%PLOT2D_SCATTER_EFFM  Similar to plot2d_scatter but highlights eff. dens. mode. 
% Author: Timothy Sipkens, 2020-11-16
%=========================================================================%

function [h] = plot2d_scatter_effm(b, grid_b, pha_rho, grid_rho)

% Get local copy of elements.
elements = grid_b.elements;


% Plot data, b.
rho_b = 6 .* elements(:,1) ./ ...
    (pi .* elements(:,2).^3) .* 1e9;
sz_b = 90 .* log10(b ./ min(b + 1 + 1e-1) + 1 + 1e-1) - 3;  % size of data pts.
scatter(elements(:,2), rho_b, ...
    sz_b, [0.85,0.85,0.85], 'filled');  % generate scatter plot


% Highlight modes of mass-selected data (e.g., mobility scan).
[~, max_b] = max(grid_b.reshape(b), [], 2);
t1=[]; t2=[]; t0=grid_b.edges{2}; sz_b_rs = grid_b.reshape(sz_b);
for ii=1:length(max_b); t1(ii)=t0(max_b(ii)); t2(ii)=sz_b_rs(ii,max_b(ii)); end
hold on;
h = scatter(t1, ...
    6 .* grid_b.edges{1} ./ ...
    (pi .* t1.^3) .* 1e9, ...
    t2, [1,0,0]);  % plot as red, empty circle
hold off;

% set axis to log-log. 
set(gca, 'XScale', 'log', 'YScale', 'log'); 


% Overlay phantom if relevant.
if exist('grid_rho', 'var')
    tools.overlay_phantom(pha_rho, grid_rho, 2, 'k');
end


% Suppress output if not requested.
if nargout==0; clear h; end

end

