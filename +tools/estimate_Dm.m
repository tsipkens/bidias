
% ESTIMATE_DM  Estimate the mass-mobility exponent from the data peaks.
% Author: Timothy Sipkens, 2020-07-17
%=========================================================================%

function [Dm] = estimate_Dm(grid_b, b)

if isempty(grid_b); Dm = []; return; end


b_rs = grid_b.reshape(b);


% fit a line to the peaks
[max_b, idx_max] = max(b_rs');
% L = diag(1./sqrt(max_b)); % weight data
A = [log10(grid_b.edges{2}(idx_max))',...
    ones(size(idx_max))']; % linear regression matrix
x = ((A)'*(A)) \ (A'*(log10(grid_b.edges{1}')));


Dm = x(1);

end

