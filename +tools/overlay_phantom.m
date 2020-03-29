
% OVERLAY_PHANTOM  Overlay circles for std. dev. isolines for phantom (optionally, fits a unimodal phantom).
% Author: Timothy Sipkens, 2019-10-31
%=========================================================================%

function [pha,N] = overlay_phantom(x_pha,grid,iso_levels,varargin)

%-- Parse inputs ---------------------------------------------------------%
if isa(x_pha,'Phantom')
    pha = x_pha;
    grid = pha.grid;
    N = [];
else
    [pha,N] = Phantom.fit(x_pha,grid);
end

if ~exist('iso_levels','var'); iso_levels = []; end
if isempty(iso_levels); iso_levels = [1,2,3]; end % plot 1, 2, and 3 sigma

if isempty(varargin); varargin = {'Color',[1,1,0,1]}; end
    % specify line properties (default: yellow, no transparency)
%-------------------------------------------------------------------------%


ylim_st = ylim;
xlim_st = xlim;

%-- Proceed with plotting ----------------%
mu = pha.mu; % local copies of parameters
Sigma = pha.Sigma;
p = pha.p;
if ~iscell(mu); mu = {mu}; Sigma = {Sigma}; end % handle unimodal case
for ii=1:pha.n_modes
    for jj=1:length(iso_levels)
        tools.overlay_ellipse(mu{ii},Sigma{ii},iso_levels(jj),varargin{:});
    end
    tools.overlay_line(grid,fliplr(mu{ii}),p(ii).Dm,varargin{:});
end


ylim(ylim_st); % restore original plot bounds
xlim(xlim_st);

if nargout==0; clear pha N; end % prevent unnecessary output

end

