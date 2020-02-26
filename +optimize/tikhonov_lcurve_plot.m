
% TIKHONOV_LCURVE_PLOT  Plot the L-curve for Tikhonov data.
% Author: Timothy Sipkens, 2020-02-06
% 
% Note: 
%   Requires 'output' structure from tikhonov_op function and plots
%   on lambda vector specified there.
%=========================================================================%

function [] = tikhonov_lcurve_plot(A,b,output)

LAxb = zeros(size(output));
x_norm = zeros(size(output));
for ii=1:length(output) % loop through pre-processed Tikhonov results
    LAxb(ii) = norm(A*output(ii).x-b)^2;
    x_norm(ii) = norm(output(ii).x);
end

loglog(LAxb,x_norm,'o');
text(LAxb,x_norm,num2cell([output.lambda]),...
    'VerticalAlignment','bottom','HorizontalAlignment','right');

end
