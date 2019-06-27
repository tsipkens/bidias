
% PLOT2D    Plots x on the grid.
% Author:   Timothy Sipkens, 2018-11-21
%=========================================================================%

function [h,x] = plot2d(obj,x)

x = reshape(x,obj.ne);

if strcmp('linear',obj.discrete)
    h = imagesc(obj.edges{2},obj.edges{1},x);
    set(gca,'YDir','normal');
    
elseif strcmp('logarithmic',obj.discrete)
    h = imagesc(log10(obj.edges{2}),log10(obj.edges{1}),x);
    set(gca,'YDir','normal');
end

end

