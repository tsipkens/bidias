
function [Lambda] = tfer_CPMA_tri(m_star,m,d,z,prop,varargin)
% TFER_CPMA_A Evaluates the transfer function for a CPMA as a triangular function.
% Author:       Timothy Sipkens, 2018-12-27
% 
%--------------------------------------------------------------------------%
% Inputs:
%   prop        CPMA properties
%   lam         CPMA geometry parameter
%   rs          Equilibirum radius
%
% Outputs:
%   Lambda      CPMA transfer function
%   G0          Function mapping final to initiral radial position
%-------------------------------------------------------------------------%


kernel_CPMA.get_setpoint; % get setpoint

if ~isfield(sp,'m_max') % if m_max/resolution was not specified
    n_B = -0.6436;
    B_star = kernel_CPMA.mp2zp(m_star,1,prop.T,prop.p); % involves invoking mass-mobility relation
    
    omega = sp.omega1.*...
        ((prop.r_hat^2-prop.omega_hat)/(prop.r_hat^2-1)+...
        prop.r1^2*(prop.omega_hat-1)/(prop.r_hat^2-1)/prop.rc^2);
    
    m_ratio = fzero(@(x) abs(x).^(n_B+1)-abs(x).^n_B-...
        prop.Q/(omega.^2.*m_star*B_star*2*pi*prop.rc^2*prop.L),m_star.*1.5);
    sp.m_max = m_star.*abs(m_ratio);
    
end


m_del = sp.m_max-m_star;
m_min = 2.*m_star-sp.m_max;

Lambda = zeros(size(m))+...
    (m<=m_star).*(m>m_min).*(m-m_min)./m_del+...
    (m>m_star).*(m<sp.m_max).*((m_star-m)./m_del+1);

end

