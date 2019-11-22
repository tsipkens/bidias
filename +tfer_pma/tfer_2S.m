
% TFER_2S   Evaluates the transfer function for a PMA in Case C.
% Author:   Timothy Sipkens, 2019-03-21
%-------------------------------------------------------------------------%
% Inputs:
%   sp          Structure defining various setpoint parameters 
%               (e.g. m_star, V). Use 'get_setpoint' method to generate 
%               this structure.
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length)
%
% Outputs:
%   Lambda      Transfer function
%   G0          Function mapping final to initial radial position
%=========================================================================%

function [Lambda,G0] = tfer_2S(sp,m,d,z,prop)

[tau,~,~,rs] = tfer_pma.parse_inputs(sp,m,d,z,prop);
        % parse inputs for common parameters

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

