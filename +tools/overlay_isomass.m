
% ISOMASS Overlay mass isolines on effective density-mobility plots.
% Author: Timothy Sipkens, 2019-11-01
%=========================================================================%

function [] = overlay_isomass(grid,color)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('color','var'); color = []; end
if isempty(color)
    alpha = 1;
    color = [1,1,1,alpha];
end
%-------------------------------------------------------------------------%


isom_vec = 10.^(-4:1:3);

for ii=1:length(isom_vec)
    hl = grid.overlay_line(...
            [1,log10(6*isom_vec(ii)/(pi*10^3).*1e9)],-3); % line for rho_eff = 0.1
    hl.Color = color;
end

end