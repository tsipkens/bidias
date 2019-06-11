
function [Lambda,G0] = tfer_CPMA_C_diff(m_star,m,d,z,prop,varargin)
% TFER_CPMA_C_DIFF Evaluates the transfer function for a CPMA in Case C.
% Author:       Timothy Sipkens, 2018-12-27
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

if ~exist('d','var')
    B = kernel_CPMA.mp2zp(m,z,prop.T,prop.p);
else
    B = kernel_CPMA.dm2zp(d,z,prop.T,prop.p);
end

D = prop.D(B).*z;
sig = sqrt(2.*prop.L.*D./prop.v_bar);

[~,G0] = kernel_CPMA.tfer_CPMA_C(m_star,m,d,z,prop,varargin{:});

rho_fun = @(G,r) (G-r)./(sqrt(2).*sig);
K_fun = @(G,r) ...
    (G-r).*erf(rho_fun(G,r))+...
    sig.*sqrt(2/pi).*exp(-rho_fun(G,r).^2);

K22 = K_fun(G0(prop.r2),prop.r2);
K21 = K_fun(G0(prop.r2),prop.r1);
K12 = K_fun(G0(prop.r1),prop.r2);
K11 = K_fun(G0(prop.r1),prop.r1);
Lambda = max(-1/(4*prop.del).*(K22-K12-K21+K11),0);


end

