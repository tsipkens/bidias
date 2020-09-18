
% OVERLAY_ELLIPSE  Plots an ellipse given a 2D mean and covariance.
% Author: Timothy Sipkens, 2019-10-27
% 
%-------------------------------------------------------------------------%
% Inputs:
%   mu        Mean as a 1x2 vector
%   Sigma     Covariance matrix as a 2x2 matrix
%   s         Standard deviation of ellipse to plot
%             (e.g., s = 2 plots an isoline ellipse two standard deviations
%             from the mean)
%   varargin  The standard extra formatting parameters to be passed to plot
%             (Optional, default: {'w'})
%-------------------------------------------------------------------------%
% 
%=========================================================================%

function h = overlay_ellipse(mu, Sigma, s, varargin)

if isempty(varargin); varargin = {'w'}; end
    % specify line properties

if ~exist('s','var'); s = []; end
if isempty(s); s = 1; end


mu = fliplr(mu(:)'); % flip mean and covar., (:)' ensures row format
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
