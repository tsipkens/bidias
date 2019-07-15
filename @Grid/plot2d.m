
% PLOT2D    Plots x on the grid.
% Author:   Timothy Sipkens, 2018-11-21
%=========================================================================%

function [h,x] = plot2d(obj,x)

x = reshape(x,obj.ne);

if strcmp('linear',obj.discrete)
    imagesc(obj.edges{2},obj.edges{1},x);
    set(gca,'YDir','normal');
    
    xlim(obj.span(2,:));
    ylim(obj.span(1,:));
    
elseif strcmp('logarithmic',obj.discrete)
    imagesc(log10(obj.edges{2}),log10(obj.edges{1}),x);
    set(gca,'YDir','normal');
    
    xlim(log10(obj.span(2,:)));
    ylim(log10(obj.span(1,:)));
end

if nargout>0; h = gca; end

end

