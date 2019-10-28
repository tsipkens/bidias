
% TFER_2S   Evaluates the transfer function for a PMA in Case C.
% Author:   Timothy Sipkens, 2019-03-21
%=========================================================================%

function [Lambda,G0] = tfer_2S(m_star,m,d,z,prop,varargin)
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


[sp,tau,C0] = ...
    tfer_pma.get_setpoint(m_star,m,d,z,prop,varargin{:});
        % get setpoint (parses d and z)

%-- Estimate equilibrium radius ------------------------------------------%
if round((sqrt(C0./m_star)-sqrt(C0./m_star-4*sp.alpha*sp.beta))/(2*sp.alpha),15)==prop.rc
    rs = real((sqrt(C0./m)-sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
else
    rs = real((sqrt(C0./m)+sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
end


%-- Estimate device parameter --------------------------------------------%
lam = 2.*tau.*(sp.alpha^2-sp.beta^2./(rs.^4)).*prop.L./prop.v_bar;


%-- Taylor series expansion constants ------------------------------------%
C1 = 2.*tau.*(sp.alpha^2-sp.beta^2./(rs.^4));
C2 = -2.*tau.*(sp.alpha^2./rs-5*sp.beta^2./(rs.^5));


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) rs+C1.*(r-rs).*exp(-lam)./...
    (C2.*(r-rs)+C1-C2.*(r-rs).*exp(-lam));

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = (1/(2*prop.del)).*(rb-ra);

end

