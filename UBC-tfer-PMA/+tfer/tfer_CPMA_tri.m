
function [Lambda] = tfer_CPMA_tri(m_star,m,d,z,prop,varargin)
% TFER_CPMA_TRI Evaluates the transfer function for a CPMA as a triangular function.
% Author: Timothy Sipkens, 2018-12-27
% 
%-------------------------------------------------------------------------%
% Inputs:
%   m_star      Setpoint particle mass
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Properties of particle mass analyzer
%   varargin    Name-value pairs for setpoint    (Optional, default Rm = 3)
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Output:
%   Lambda      Transfer function
%-------------------------------------------------------------------------%


tfer.get_setpoint; % get setpoint

if ~isfield(sp,'m_max') % if m_max was not specified
    n_B = -0.6436;
    B_star = tfer.mp2zp(m_star,1,prop.T,prop.p); % involves invoking mass-mobility relation
    
    omega = sp.omega1.*...
        ((prop.r_hat^2-prop.omega_hat)/(prop.r_hat^2-1)+...
        prop.r1^2*(prop.omega_hat-1)/(prop.r_hat^2-1)/prop.rc^2);
    
    m_ratio = fzero(@(x) abs(x).^(n_B+1)-abs(x).^n_B-...
        prop.Q/(omega.^2.*m_star*B_star*2*pi*prop.rc^2*prop.L),m_star.*1.5);
    sp.m_max = m_star.*abs(m_ratio);
    
end


m_del = sp.m_max-m_star; % FWHM of the transfer function (related to resolution)
m_min = 2.*m_star-sp.m_max; % lower end of the transfer function

Lambda = zeros(size(m))+...
    (m<=m_star).*(m>m_min).*(m-m_min)./m_del+...
    (m>m_star).*(m<sp.m_max).*((m_star-m)./m_del+1);
        % evaluate the transfer function

end

