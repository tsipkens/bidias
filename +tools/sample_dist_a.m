
% SAMPLE_DIST_A  Sample from a discrete 2D distribution using 'slicesample'.
% Author: Timothy Sipkens, 2020-03-30
%=========================================================================%

function [rnd] = sample_dist_a(x,grid,nsamples)

if~exist('nsamples','var'); nsamples = []; end
if isempty(nsamples); nsamples = 1e4; end

[~,ind_max] = max(x);
initial = log10(grid.elements(ind_max,:));

rnd = slicesample(initial,nsamples,'pdf',@(r) pdf_fun(x,grid,r));

end


function p = pdf_fun(x,grid,r)

k = grid.closest_idx(10.^r);

p = zeros(size(r,1),1);
p(~isnan(k)) = x(k(~isnan(k)));

end