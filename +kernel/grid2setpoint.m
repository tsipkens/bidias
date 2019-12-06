
% GRID2SETPOINT Convert a grid used for data to a structured array of setpoints.
% Author: Timothy Sipkens, 2019-12-06
%=========================================================================%

function [sp] = grid2setpoint(prop,grid_b,Rm)

m_star = grid_b.elements(:,1);

for ii=length(m_star):-1:1
    sp(ii) = tfer_pma.get_setpoint(prop,...
        'm_star',m_star(ii),'Rm',Rm);
end

end

