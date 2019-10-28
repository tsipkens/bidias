
% TFER_GE   Evaluates the transfer function for a PMA in Case F.
% Author:   Timothy Sipkens, 2019-03-21
%=========================================================================%

function [Lambda,G0] = tfer_GE(m_star,m,d,z,prop,varargin)
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


%-- Estimate recurring quantities ----------------------------------------%
C6 = 2*sp.alpha*sp.beta-C0./m;
C7 = sqrt(4*sp.alpha^2*sp.beta^2-C6.^2);

A1 = prop.v_bar./(4.*tau.*sp.alpha^2);
A2 = 2.*C6./C7;


%-- Set up F function for minimization -----------------------------------%
F = @(r,ii) A1(ii).*(log(sp.alpha^2.*r^4+sp.beta^2+C6(ii).*r.^2)-...
    A2(ii).*atan((2*sp.alpha^2.*r.^2+C6(ii))./C7(ii)));
min_fun = @(rL,r0,ii) F(rL,ii)-F(r0,ii)-prop.L;


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) tfer_pma.G_fun(min_fun,r,rs,prop.r1,prop.r2,sp.alpha,sp.beta);

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = (1/(2*prop.del)).*(rb-ra);

end

