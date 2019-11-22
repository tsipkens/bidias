
% TFER_TRI Evaluates the transfer function for a PMA as a triangular function.
% Author: Timothy Sipkens, 2018-12-27
%%-------------------------------------------------------------------------%
% Inputs:
%   sp          Structure defining various setpoint parameters 
%               (e.g. m_star, V). Use 'get_setpoint' method to generate 
%               this structure.
%   m           Particle mass
%   ~           Placeholder for mobility diameter
%   ~           Placeholder for integer charge state
%   ~           Placeholder for device properties (e.g. classifier length)
%
% Outputs:
%   Lambda      Transfer function
%=========================================================================%

function [Lambda] = tfer_tri(sp,m,~,~,~)

m_del = sp.m_max-sp.m_star; % FWHM of the transfer function (related to resolution)
m_min = 2.*sp.m_star-sp.m_max; % lower end of the transfer function

Lambda = zeros(size(m))+...
    (m<=sp.m_star).*(m>m_min).*(m-m_min)./m_del+...
    (m>sp.m_star).*(m<sp.m_max).*((sp.m_star-m)./m_del+1);
        % evaluate the transfer function

end

