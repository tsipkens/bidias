
function [Lambda,G0] = tfer_CPMA_E_pb(m_star,m,d,z,prop,varargin)
% TFER_CPMA_E_PB Evaluates the transfer function for a CPMA in Case E.
% Author:       Timothy Sipkens, 2019-03-21
% 
%-------------------------------------------------------------------------%
% Inputs:
%   m_star      Setpoint particle mass
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Properties of the particle parameters
%   varargin    Name-value pairs for setpoint    (Optional, default Rm = 3)
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Outputs:
%   Lambda      CPMA transfer function
%   G0          Function mapping final to initiral radial position
%-------------------------------------------------------------------------%


kernel_CPMA.get_setpoint; % get setpoint

%-- Estimate equilibrium radius ------------------------------------------%
if round((sqrt(C0./m_star)-sqrt(C0./m_star-4*sp.alpha*sp.beta))/(2*sp.alpha),15)==prop.rc
    rs = real((sqrt(C0./m)-sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
else
    rs = real((sqrt(C0./m)+sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
end


%-- Estimate recurring quantities ----------------------------------------%
A1 = -3.*prop.v_bar./(4.*tau.*sp.omega1.^4.*prop.del^2);
A2 = sp.omega1.^2.*(prop.rc^2-prop.del^2)+C0./m;
A3 = 4*prop.rc.*sp.omega1.*sqrt(C0./m);
A4 = @(r,ii) sp.omega1.^2.*(r.^2-4*prop.rc.*r);


%-- Set up F function for minimization -----------------------------------%
F = @(r,ii) A1(ii).*(A2(ii).*log(C0./m(ii)-(sp.omega1.*r).^2)+...
    A3(ii).*atanh(sqrt(m(ii)./C0).*sp.omega1.*r)+A4(r,ii));
min_fun = @(rL,r0,ii) F(rL,ii)-F(r0,ii)-prop.L;

condit = (sp.alpha^2) < (sp.beta^2./(rs.^4));


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) kernel_CPMA.G_fun(min_fun,r,rs,prop.r1,prop.r2,condit);

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = 3/4.*(rb-ra)./prop.del-1/4.*((rb-prop.rc)./prop.del).^3+...
    1/4.*((ra-prop.rc)./prop.del).^3;

end

