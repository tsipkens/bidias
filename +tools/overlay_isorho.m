
% OVERLAY_ISORHO Overlay effective density isolines on mass-mobility plots.
% Author: Timothy Sipkens, 2019-11-01
%=========================================================================%

function [] = overlay_isorho(grid,color)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('color','var'); color = []; end
if isempty(color)
    alpha = 1;
    color = [1,1,1,alpha];
end
%-------------------------------------------------------------------------%


isorho_vec = 10.^(-1:1:7);
for ii=1:length(isorho_vec)
    hl = grid.overlay_line(...
            [1,log10(isorho_vec(ii)*10^3./1e9)],3); % line for rho_eff = 0.1
    hl.Color = color;
end

end