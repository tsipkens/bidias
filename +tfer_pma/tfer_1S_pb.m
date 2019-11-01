
% TFER_1S_PB    Evaluates the transfer function for a PMA in Case A (w/ parabolic flow).
% Author:       Timothy Sipkens, 2018-12-27
%=========================================================================%

function [Lambda,G0,sp] = tfer_1S_pb(m_star,m,d,z,prop,varargin)
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

A1 = -3*prop.L./(2.*lam.*prop.del^2);
A2 = rs.^2+prop.rc^2-2*prop.rc.*rs-prop.del^2;
A3 = @(r,ii) (r.^2)./2+(rs(ii)-2*prop.rc).*r;


%-- Set up F function for minimization -----------------------------------%
F = @(r,ii) A1(ii).*(A2(ii).*log(r-rs(ii))+A3(r,ii));
min_fun = @(rL,r0,ii) F(rL,ii)-F(r0,ii)-prop.L;


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) tfer_pma.G_fun(min_fun,r,rs,prop.r1,prop.r2,sp.alpha,sp.beta);

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = 3/4.*(rb-ra)./prop.del-1/4.*((rb-prop.rc)./prop.del).^3+...
    1/4.*((ra-prop.rc)./prop.del).^3;

end

