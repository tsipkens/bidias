
% TFER_2C_DIFF  Evaluates the transfer function for a PMA in Case D (w/ diffusion).
% Author:       Timothy Sipkens, 2018-12-27
%=========================================================================%

function [Lambda,G0] = tfer_2C_diff(m_star,m,d,z,prop,varargin)
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


%-- Evaluate mechanical mobility for diffusion calc. ---------------------%
if ~exist('d','var')
    B = tfer_PMA.mp2zp(m,z,prop.T,prop.p);
        % if mobility is not specified, use mass-mobility relation to estimate
else
    B = tfer_PMA.dm2zp(d,z,prop.T,prop.p);
end

D = prop.D(B).*z;
    % diffusion coefficient is previously defined function multiplied by 
    % integer charge state
sig = sqrt(2.*prop.L.*D./prop.v_bar); % diffusive spreading parameter

[~,G0] = tfer_PMA.tfer_2C(m_star,m,d,z,prop,varargin{:});
    % get G0 function for this case

rho_fun = @(G,r) (G-r)./(sqrt(2).*sig); % reuccring quantity
kap_fun = @(G,r) ...
    (G-r).*erf(rho_fun(G,r))+...
    sig.*sqrt(2/pi).*exp(-rho_fun(G,r).^2); % define function for kappa

%-- Evaluate the transfer function and its terms -------------------------%
K22 = kap_fun(real(G0(prop.r2)),prop.r2);
K21 = kap_fun(real(G0(prop.r2)),prop.r1);
K12 = kap_fun(real(G0(prop.r1)),prop.r2);
K11 = kap_fun(real(G0(prop.r1)),prop.r1);
Lambda = -1/(4*prop.del).*(K22-K12-K21+K11);
Lambda(K22>1e2) = 0; % remove cases with large values out of error fun. eval.
Lambda(abs(Lambda)<1e-10) = 0; % remove cases with roundoff error

end
