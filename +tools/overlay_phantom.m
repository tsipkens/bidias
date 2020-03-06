
% OVERLAY_PHANTOM Fit a distribution and overlay circles for std. dev. isolines
% Author: Timothy Sipkens, 2019-10-31
%=========================================================================%

function [pha,N] = overlay_phantom(x_pha,grid,varargin)

ylim_st = ylim;
xlim_st = xlim;

if isa(x_pha,'Phantom')
    pha = x_pha;
    grid = pha.grid;
    N = [];
else
    [pha,N] = Phantom.fit(x_pha,grid);
end

if isempty(varargin); varargin = {'Color',[1,1,0,1]}; end
    % specify line properties (default: yellow, no transparency)

%-- Proceed with plotting ----------------%
tools.overlay_ellipse(pha.mu,pha.Sigma,1,varargin{:});
tools.overlay_ellipse(pha.mu,pha.Sigma,2,varargin{:});
tools.overlay_ellipse(pha.mu,pha.Sigma,3,varargin{:});

tools.overlay_line(grid,fliplr(pha.mu),...
    pha.p.Dm,varargin{:});

ylim(ylim_st); % restore original plot bounds
xlim(xlim_st);

if nargout==0; clear pha N; end % prevent unnecessary output

end

