
% OVERLAY_ELLIPSE Plots an ellipse give a center, major and minor radii, and slope.
% Author: Timothy Sipkens, 2019-10-27
% Note:   For raidus, the first entry is the major axis, the second entry
%         is the minor axis
%=========================================================================%

function h = overlay_ellipse(mu,Sigma,s,cspec)

if ~exist('cspec','var'); cspec = []; end
if isempty(cspec); cspec = 'w'; end

mu = fliplr(mu); % flip mean and covar.
Sigma = rot90(Sigma,2);

s = s.*2; % double number of std. dev.

[V, D] = eig(Sigma.*s);

t = linspace(0, 2 * pi);
a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

hold on;
h = plot(a(1, :) + mu(1), a(2, :) + mu(2), cspec);
hold off;


%{
s = -2 * log(1 - p);
[V,D] = eig(Sigma*s);

theta = 0:0.01:2*pi;
y = Sigma(2,2) * cos(theta);
x = Sigma(1,1) * sin(theta);
xy = [x;y]';


%-- Transform circle based on slope --------------------------------------%
rotate = atan(sigma(1)/sigma(2))*180/pi;
transform_matrix = [cosd(rotate), sind(rotate);...
  -sind(rotate), cosd(rotate)];
xy_rotate = xy * transform_matrix;
x_rotate = xy_rotate(:, 1) + center(2);
y_rotate = xy_rotate(:, 2) + center(1);


%-- Plot the ellipse -----------------------------------------------------%
hold on;
plot(center(2),center(1),[cspec,'x']); % plot center of ellipse
h = plot(x_rotate,y_rotate,cspec); % plot ellipse
hold off;

%}

if nargout==0; clear h; end % remove output if none requested

end
