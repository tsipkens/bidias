
% OVERLAY_PEAKS  Use a traditional analysisto plto peaks on 2D distribution
%  
%  AUTHOR: Timothy Sipkens, 2022-10-20

function [] = overlay_peaks(grid, x, dim, cspec)

dim2 = setdiff([1,2],dim);

x_rs = grid.reshape(x);

[~, idx] = max(x_rs, [], dim2);

% Filter out edges, when found.
fidx = true(size(idx));
fidx(idx == 1) = 0;
fidx(idx == size(x_rs, dim)) = 0;
idx(idx == 1) = [];
idx(idx == size(x_rs, dim)) = [];

if dim == 2
    plot(grid.edges{2}(fidx), grid.edges{1}(idx), cspec);
else
    plot(grid.edges{2}(idx), grid.edges{1}(fidx), cspec);
end

end

