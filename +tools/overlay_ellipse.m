
% OVERLAY_ELLIPSE Plots an ellipse give a center, major and minor radii, and slope.
% Author: Timothy Sipkens, 2019-10-27
% Note:   For raidus, the first entry is the major axis, the second entry
%         is the minor axis
%=========================================================================%

function h = overlay_ellipse(mu,Sigma,s,varargin)

if isempty(varargin); varargin = {'w'}; end
    % specify line properties

if ~exist('s','var'); s = []; end
if isempty(s); s = 1; end


mu = fliplr(mu); % flip mean and covar.
Sigma = rot90(Sigma,2);

s = s.*2; % double number of std. dev. for ellipse
          % (i.e. one std. dev. in each direction)

[V, D] = eig(Sigma.*s);

t = linspace(0, 2*pi);
a = (V*sqrt(D))*[cos(t(:))'; sin(t(:))'];

hold on;
if strcmp(get(gca,'XScale'),'log')
    h = loglog(10.^(a(1,:)+mu(1)), 10.^(a(2,:)+mu(2)), varargin{:});
else
    h = plot(a(1,:)+mu(1), a(2,:)+mu(2), varargin{:});
end
hold off;

if nargout==0; clear h; end % remove output if none requested

end
