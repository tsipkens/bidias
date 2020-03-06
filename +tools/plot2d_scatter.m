
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
if isempty(cmap); cmap = colormap('gray'); end

N = size(cmap,1);
color = cmap(round(corder.*(N-1)+1),:);
color(corder==1,:) = 1;

clf;
for ii=1:length(m)
    if ii==2; hold on; end
    loglog(d(ii),m(ii),'.',...
        'Color',color(ii,:),...
        'MarkerSize',marker_size(ii));
        % marker_size is also logscale
end
hold off;

%-- Add custom legend ----------------------%
hold on;
for ii=5:-1:0
    t0 = (bmax+ii-5-bmin)/(bmax-bmin);
    if t0>=0
        h(6-ii) = plot(NaN,NaN,'.',...
            'Color',cmap(round((1-t0)*(N-1)+1),:),...
            'MarkerSize',12*t0+0.1);
        text{6-ii} = num2str(10^(ii-5));
    end
end
text{end} = ['<',text{end}];
hold off;

lgd = legend(h,text,'Location','northwest');
title(lgd,'{\times} b_{max}');
%-------------------------------------------%

end