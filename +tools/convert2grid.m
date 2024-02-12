
% CONVERT2GRID  A utility to convert a list of setpoints to a grid.
% 
% Author: Timothy Sipkens, 2020-07-16
%=========================================================================%

function [grid, idx, sp] = convert2grid(sp_m, d_star, b)

% adjust mass setpoint
if isstruct(sp_m); m_star = [sp_m.m_star] .* 1e18; % if setpoint structure
else; m_star = sp_m; end % if mass setpoints

b_size = length(d_star);


% round the setpoints to discretize
m_star = exp(round(log(m_star), 1));  % grouped based on one-digit rounding
d_star = exp(round(log(d_star), 2))'; % grouped based on two-digit rounding

[m_star, ia_m, ic_m] = unique(m_star);
[d_star, ~, ic_d] = unique(d_star);
[~,idx] = sortrows([ic_d, ic_m]);

if length(m_star) * length(d_star) == b_size
    grid = Grid({m_star, d_star});
    sp = sp_m(ia_m);

else
    disp('Failed to convert to Grid.');
    grid = [];
    idx = [];
end

end

