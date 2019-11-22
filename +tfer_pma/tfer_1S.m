
% TFER_1S   Evaluates the transfer function for a PMA in Case A.
% Author:   Timothy Sipkens, 2018-12-27
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

function [Lambda,G0] = tfer_1S(sp,m,d,z,prop)

[tau,~,~,rs] = tfer_pma.parse_inputs(sp,m,d,z,prop);
        % parse inputs for common parameters

%-- Estimate device parameter --------------------------------------------%
lam = 2.*tau.*(sp.alpha^2-sp.beta^2./(rs.^4)).*prop.L./prop.v_bar;


%-- Evaluate G0 and transfer function ------------------------------------%
G0 = @(r) rs+(r-rs).*exp(-lam);

ra = min(prop.r2,max(prop.r1,G0(prop.r1)));
rb = min(prop.r2,max(prop.r1,G0(prop.r2)));

Lambda = (1/(2*prop.del)).*(rb-ra);

end

