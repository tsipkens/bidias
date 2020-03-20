
% OVERLAY_PHANTOM_MARG  Overlay ellipses and plot on marginal panels (optionally, fits a unimodal phantom). 
% Author: Timothy Sipkens, 2020-03-20
%=========================================================================%

function [pha,N] = overlay_phantom_marg(x_pha,N,grid,iso_levels,varargin)

% Parse inputs: Perform initial ellipse overlays for phantom + optionally fit unimodal phantom
if isa(x_pha,'Phantom') % don't fit phantom
    pha = x_pha;
    tools.overlay_phantom(pha,grid,iso_levels,varargin{:});
    
else % fit phantom
    [pha,N] = tools.overlay_phantom(x_pha,grid,iso_levels,varargin{:});
end


%-- Generate higher resolution plot vectors ------------------------------%
subplot(4,4,[1,3]); % marginal panel for dim2
xlim_st = xlim;
vec_dim2 = logspace(log10(xlim_st(1)),log10(xlim_st(2)),200);

subplot(4,4,[8,16]); % marginal panel for dim1
ylim_st = ylim;
vec_dim1 = logspace(log10(ylim_st(1)),log10(ylim_st(2)),199);

pha_dim1 = zeros(size(vec_dim1));
pha_dim2 = zeros(size(vec_dim2));


%-- Evaluate phantom marginals -------------------------------------------%
for ii=1:pha.n_modes
    pha_dim1 = pha_dim1 + ...
        N(ii).*normpdf(log10(vec_dim1),...
        log10(pha.p(ii).mg),log10(pha.p(ii).sm));
    
    pha_dim2 = pha_dim2 + ...
        N(ii).*normpdf(log10(vec_dim2),...
        log10(pha.p(ii).dg),log10(pha.p(ii).sg));
end


%-- Plot on marginal panels ----------------------------------------------%
subplot(4,4,[1,3]); % add phantom to mobility (or dim2) marginal dist.
hold on;
semilogx(vec_dim2,pha_dim2);
hold off;

subplot(4,4,[8,16]); % add phantom to mass (or dim1) marginal dist.
hold on;
semilogx(pha_dim1,vec_dim1);
hold off;



subplot(4,4,[5,15]); % refocus center panel

if nargout==0; clear pha N; end % prevent unnecessary output

%-------------------------%
%{
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
%}

end

