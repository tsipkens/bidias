
% REDISTRIBUTE  Estimate the distribution width of a transformed random variable.
% Author: Timothy Sipkens, 2020-03-21
%=========================================================================%

function s = redistribute(s0,fun,mu)

if ~exist('mu','var'); mu = []; end
if isempty(mu); mu = 0; end

n = 6000;
s = std(fun(mu+s0.*randn([size(s0),n])),[],ndims(s0)+1);

end

