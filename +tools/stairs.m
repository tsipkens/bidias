
% STAIRS  A modifed stairs function for histograms. 
%  
%  AUTHOR: Timothy Sipkens, 2022-10-20

function p = stairs(x, y, flog, varargin)

if ~exist('flog', 'var'); flog = []; end
if isempty(flog); flog = 0; end

if flog
    mx = exp(log(x(2:end)) + log(x(1:(end-1))) ./ 2);  % midpoints in x
    dx = exp(diff(log(x)));
else
    mx = (x(2:end) + x(1:(end-1))) ./ 2;  % midpoints in x
    dx = diff(x);
end

p = stairs([x(1) - dx(1), mx, x(end) + dx(end)], [y, 0], varargin{:});

if nargout == 0; clear p; end

end