
% SAMPLE_DIST_c  Sample from a discrete 2D distribution using 'slicesample'.
% Author: Timothy Sipkens, 2020-03-30
%=========================================================================%

function [rnd] = sample_dist_c(x,vec,nsamples)

if~exist('nsamples','var'); nsamples = []; end
if isempty(nsamples); nsamples = 1e4; end

[~,ind_max] = max(x);
initial = log10(vec(ind_max,:));
F = scatteredInterpolant(log10(vec(:,1)),log10(vec(:,2)),x,...
    'linear','none');

rnd = slicesample(initial,nsamples,'pdf',@(r) pdf_fun(F,r));

end


function p = pdf_fun(F,r)

p = F(r(:,1),r(:,2));
p(isnan(p)) = 0;

end