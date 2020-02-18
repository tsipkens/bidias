
% TFER_EHARA    Evaluates the transfer function for a PMA in Case 1S as per Ehara et al. (1996).
% Author:       Timothy Sipkens, 2018-12-27
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
%=========================================================================%

function [Lambda] = tfer_ehara(sp,m,d,z,prop)

[tau,~,~,rs] = tfer_pma.parse_inputs(sp,m,d,z,prop);
        % parse inputs for common parameters

%-- Estimate device parameter --------------------------------------------%
lam = 2.*tau.*(sp.alpha^2-sp.beta^2./(rs.^4)).*prop.L./prop.v_bar;


%-- Evaluate transfer function -------------------------------------------%
rho_s = (rs-prop.rc)/prop.del;
Lambda = ((1-rho_s)+(1+rho_s).*exp(-lam))./2.*and(1<rho_s,rho_s<coth(lam./2))+...
    exp(-lam).*and(-1<rho_s,rho_s<1)+...
    ((1+rho_s)+(1-rho_s).*exp(-lam))./2.*and(-coth(lam./2)<rho_s,rho_s<-1);


end

