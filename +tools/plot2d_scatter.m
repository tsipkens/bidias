
% PLOT2D_SCATTER Scatter plot of data, where colour indicates counts.
% Author: Timothy Sipkens, 2019-11-28
%=========================================================================%

function [] = plot2d_scatter(m,d,b,cmap)

b = b./max(b);

logb = log10(b); % scale data
bmax = max(logb);
bmin = min(logb(~isinf(logb)));
bmin = max(bmin,bmax-5); % at most, span four orders of magnitude
logb = max(logb,bmin);
bscl = max((logb-bmin)/(bmax-bmin),bmin);

marker_size = 12.*bscl+0.1;
corder = 1-bscl; % color is logscale

if ~exist('cmap','var'); cmap = []; end
if isempty(cmap)
    color = corder*ones(1,3);
else
    N = size(cmap,1);
    color = cmap(round(corder.*(N-1)+1),:);
    color(corder==1,:) = 1;
end

clf;
hold on;
for ii=1:length(m)
    plot(log10(d(ii)),log10(m(ii)),'.',...
        'Color',color(ii,:),...
        'MarkerSize',marker_size(ii));
        % marker_size is also logscale
end
hold off;

end