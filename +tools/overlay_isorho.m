
% OVERLAY_ISORHO Overlay effective density isolines on mass-mobility plots.
% Author: Timothy Sipkens, 2019-11-01
%=========================================================================%

function [] = overlay_isorho(grid, varargin)

isorho_vec = 10.^(-1:1:9);  % vector of effective densities to plot
for ii=1:length(isorho_vec)
    tools.overlay_line(grid,...
    	[1, log10(pi/6*isorho_vec(ii)*10^3./1e9)], ...  % point at d = 10 nm
        3, ...  % slope fis always 3 (eff. dens. is a cubic function of diameter)
        varargin{:});  % additional arguments (e.g., colour)
end

end