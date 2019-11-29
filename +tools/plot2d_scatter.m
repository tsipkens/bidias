
% PLOT2D_SCATTER Scatter plot of data, where colour indicates counts.
% Author: Timothy Sipkens, 2019-11-28
%=========================================================================%

function [] = plot2d_scatter(m,d,b,cmap)

b = b./max(b);

corder = min(1-1/3.*(log10([b,b,b])+3),1);

clf;
hold on;
for ii=1:length(m)
    plot(log10(d(ii)),log10(m(ii)),'.',...
        'Color',corder(ii,:),...
        'MarkerSize',max((log10(b(ii))+8),1));
end
hold off;

end