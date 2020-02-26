
% TIKHONOV_LCURVE_PLOT  Plot the L-curve for Tikhonov data.
% Requires 'out' structure from tikhonov_op function.
% Author: Timothy Sipkens, 2020-02-06
%=========================================================================%

function [] = tikhonov_lcurve_plot(A,b,out)

LAxb = zeros(size(out));
x_norm = zeros(size(out));
for ii=1:length(out) % loop through pre-processed Tikhonov results
    LAxb(ii) = norm(A*out(ii).x-b)^2;
    x_norm(ii) = norm(out(ii).x);
end

loglog(LAxb,x_norm,'o');
text(LAxb,x_norm,num2cell([out.lambda]),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');

end
