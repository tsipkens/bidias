
% TFER_1C_PB    Evaluates the transfer function for a PMA in Case B (w/ parabolic flow).
% Author:       Timothy Sipkens, 2019-03-21
%=========================================================================%

function [Lambda,G0] = tfer_1C_pb(m_star,m,d,z,prop,varargin)
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


tfer_pma.get_setpoint; % get setpoint (parses d and z)

%-- Taylor series expansion constants ------------------------------------%
C3 = tau.*(sp.alpha^2*prop.rc+2*sp.alpha*sp.beta/prop.rc+sp.beta^2/(prop.rc^3)-C0./(m.*prop.rc));
C4 = tau.*(sp.alpha^2-2*sp.alpha*sp.beta/(prop.rc^2)-3*sp.beta^2/(prop.rc^4)+C0./(m.*(prop.rc^2)));

A1 = -3*prop.v_bar./(4.*C4.^3.*prop.del^2);
A2 = 2.*(C3.^2-C4.^2.*prop.del^2);
A3 = @(r,ii) C4(ii).^2.*(r-prop.rc).^2-2.*C3(ii).*C4(ii).*(r-prop.rc);


%-- Estimate equilibrium radius ------------------------------------------%
if round((sqrt(C0./m_star)-sqrt(C0./m_star-4*sp.alpha*sp.beta))/(2*sp.alpha),15)==prop.rc
    rs = real((sqrt(C0./m)-sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
else
    rs = real((sqrt(C0./m)+sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
end


%-- Set up F function for minimization -----------------------------------%
F = @(r,ii) A1(ii).*(A2(ii).*log(abs(C4(ii).*(r-prop.rc)+C3(ii)))+A3(r,ii));
min_fun = @(rL,r0,ii) F(rL,ii)-F(r0,ii)-prop.L;


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) tfer_pma.G_fun(min_fun,r,rs,prop.r1,prop.r2,sp.alpha,sp.beta);

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = 3/4.*(rb-ra)./prop.del-1/4.*((rb-prop.rc)./prop.del).^3+...
    1/4.*((ra-prop.rc)./prop.del).^3;

end

