
function [Lambda,G0] = tfer_CPMA_D(m_star,m,d,z,prop,varargin)
% TFER_CPMA_D Evaluates the transfer function for a CPMA in Case D.
% Author:       Timothy Sipkens, 2019-03-21
% 
%-------------------------------------------------------------------------%
% Inputs:
%   m_star      Setpoint particle mass
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length)
%   varargin    Name-value pairs for setpoint    (Optional, default Rm = 3)
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Outputs:
%   Lambda      Transfer function
%   G0          Function mapping final to initial radial position
%-------------------------------------------------------------------------%

kernel.get_setpoint; % get setpoint

%-- Taylor series expansion constants ------------------------------------%
C3 = tau.*(sp.alpha^2*prop.rc+2*sp.alpha*sp.beta/prop.rc+sp.beta^2/(prop.rc^3)-C0./(m.*prop.rc));
C4 = tau.*(sp.alpha^2-2*sp.alpha*sp.beta/(prop.rc^2)-3*sp.beta^2/(prop.rc^4)+C0./(m.*(prop.rc^2)));
C5 = 2.*tau.*(2*sp.alpha*sp.beta./(prop.rc^3)+6*sp.beta^2/(prop.rc^5)-C0./(m.*(prop.rc^3)));

C6 = -1./(sqrt(4.*C3.*C5-C4.^2));

f = @(r) 2.*C6.*prop.v_bar.*...
    atan(C6.*(2.*C5.*(r-prop.rc)+C4));


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) (tan((f(r)-prop.L)./(2.*C6.*prop.v_bar))./C6-C4)./(2.*C5)+prop.rc;

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = (1/(2*prop.del)).*(rb-ra);

end

