
% OVERLAY_PHANTOM Fit a distribution and overlay circles for std. dev. isolines
% Author: Timothy Sipkens, 2019-10-31
%=========================================================================%

function [pha] = overlay_phantom(x_pha,grid)

if isa(x_pha,'Phantom')
    pha = x_pha;
    grid = pha.grid;
else
    pha = Phantom.fit(x_pha,grid);
end

h1 = tools.overlay_ellipse(pha.mu,pha.Sigma,1);
h1.Color=[1,1,0,0.5];
h1 = tools.overlay_ellipse(pha.mu,pha.Sigma,2);
h1.Color=[1,1,0,0.5];
h1 = tools.overlay_ellipse(pha.mu,pha.Sigma,3);
h1.Color=[1,1,0,0.5];

h1 = grid.overlay_line(fliplr(pha.mu),pha.p.Dm);
h1.Color=[1,1,0,0.5];

if nargout==0; clear pha; end % prevent unnecessary output

end

