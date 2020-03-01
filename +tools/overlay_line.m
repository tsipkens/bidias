
% OVERLAY_LINE  Plots a line on top of the current grid
% Author:	Timothy Sipkens, 2019-07-15
%
% Inputs:
%   logr0   A single point on the line
%	slope   Slope of the line
%   c_spec  Color specification string, e.g. 'k' for a black line
% Outputs:
%   h       Line object
%=========================================================================%

function h = overlay_line(grid,logr0,slope,cspec)

if ~exist('cspec','var'); cspec = 'w'; end

if strcmp(get(gca,'XScale'),'log') % for log-scale plots
    rmin = log10(min([grid.edges{:}]));
    rmax = log10(max([grid.edges{:}]));

    hold on;
    h = loglog(10.^[rmin,rmax],...
        10.^[logr0(2)+slope*(rmin-logr0(1)),...
        logr0(2)+slope*(rmax-logr0(1))],cspec);
    hold off;

else % for linear scale plots
    rmin = min([grid.edges{:}]);
    rmax = max([grid.edges{:}]);

    hold on;
    h = plot([rmin,rmax],...
        [logr0(2)+slope*(rmin-logr0(1)),...
        logr0(2)+slope*(rmax-logr0(1))],cspec);
    hold off;
end

if nargout==0; clear h; end

end



