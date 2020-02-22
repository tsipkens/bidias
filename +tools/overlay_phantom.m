
% OVERLAY_PHANTOM Fit a distribution and overlay circles for std. dev. isolines
% Author: Timothy Sipkens, 2019-10-31
%=========================================================================%

function [pha,N] = overlay_phantom(x_pha,grid,color)

if isa(x_pha,'Phantom')
    pha = x_pha;
    grid = pha.grid;
    N = [];
else
    [pha,N] = Phantom.fit(x_pha,grid);
end

if ~exist('color','var'); color = []; end
if isempty(color); alpha = 1; color = [1,1,0,alpha]; end


h1 = tools.overlay_ellipse(pha.mu,pha.Sigma,1);
h1.Color = color;
h1 = tools.overlay_ellipse(pha.mu,pha.Sigma,2);
h1.Color = color;
h1 = tools.overlay_ellipse(pha.mu,pha.Sigma,3);
h1.Color = color;

h1 = grid.overlay_line(fliplr(pha.mu),...
    pha.p.Dm);
h1.Color = color;

if nargout==0; clear pha N; end % prevent unnecessary output

end

