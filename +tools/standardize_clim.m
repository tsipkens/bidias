
 % STANDARDIZE_CLIM  Standardize the color limit of a plot. 
 %  
 %  AUTHOR: Timothy Sipkens, 2023-03-24

function standardize_clim(clim_div)

if ~exist('f', 'var'); f = []; end
if isempty(f); f = gcf; end

if ~exist('clim_div', 'var'); clim_div = []; end

child = f.Children;

clim = [inf, -inf];
for ii=1:length(child)
    if isa(child(ii), 'matlab.graphics.axis.Axes')
        clim = [min(clim(1), child(ii).CLim(1)), ...
                max(clim(2), child(ii).CLim(2))];
    end
end

if ~isempty(clim_div)
    if length(clim_div) == 2
        clim = clim_div;

    elseif clim_div  % assume bool flagging diverging
        cmax = max(abs(clim));
        clim = [-cmax, cmax];
    end
end

for ii=1:length(child)
    if isa(child(ii), 'matlab.graphics.axis.Axes')
        set(child(ii), 'CLim', clim);
    end
end

end
