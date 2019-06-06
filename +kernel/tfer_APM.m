    
function [Lambda,prop] = tfer_APM(m_star,m,d_star,z,prop)
% F_APM     Evaluates the transfer function for a aerosol particle mass analyzer.
% Original source:  Buckley et al., J. Aerosol Sci. (2008)
% Modified by:      Timothy Sipkens, 2018-12-27
% Note: This transfer function is calculated based off Ehara et al (1996).
% 
%--------------------------------------------------------------------------%
% Inputs:
%   m_star      Mass corresponding to the measurement set point of the APM
%   d_star      Particle diameter, measurement set point for upstream DMA [nm]
%   m           Particle mass, can be vector
%   z           Integer charge state
%
% Outputs:
%   Lambda      APM transfer function
%-------------------------------------------------------------------------%

e = 1.6022E-19; % electron charge [C]


%-- APM settings ---------------------------------------------------------%
if ~exist('prop','var') % if order not specified
    prop = default_conditions;
elseif isempty(prop)
    prop = default_conditions;
end


%-- Calculate rotational speed if not specified ---------------------------%
if ~isfield(prop,'omega') % if rotational speed is not explicitly specified, Olfert lab
    rc = (prop.r2+prop.r1)/2;
    n_B = -0.6436;
    B_star = kernel.mp2zp(m_star,1,prop.T,prop.p); % use z = 1 for APM setpoint
    Rm = 10;
    m_max = m_star*(1/Rm+1);
    prop.omega = sqrt(prop.Q/(m_star*B_star*2*pi*rc^2*prop.L*...
        ((m_max/m_star)^(n_B+1)-(m_max/m_star)^n_B)));
end


%-- Transfer function calculation-----------------------------------------%
v_bar = prop.Q/(pi*(prop.r2^2-prop.r1^2)); % average axial flow through APM
rc = (prop.r1+prop.r2)/2;
del = (prop.r2-prop.r1)/2;

V = m_star.*rc^2*prop.omega^2*log(prop.r2/prop.r1)./e;
rs = (V./((m./(z*e)).*prop.omega^2*log(prop.r2/prop.r1))).^0.5;
rho = (rs-rc)/del;

if isempty(d_star)
    [~,Zp_star] = kernel.mp2zp(m_star,z,prop.T,prop.p);
else
    [~,Zp_star] = kernel.dm2zp(d_star,z,prop.T,prop.p); % Previously: kernel.invZp(d_star);
end
tau = m*Zp_star/(z*e);

lam = 2*tau*prop.omega^2*prop.L/v_bar;

Lambda = (rho>=-1).*(rho<=1).*exp(-lam)+...
    (rho>1).*(rho<coth(lam/2)).*((1-rho)+(1+rho).*exp(-lam))./2+...
    (rho<-1).*(rho>-coth(lam/2)).*((1+rho)+(1-rho).*exp(-lam))./2;
%-------------------------------------------------------------------------%

end


function prop = default_conditions()

prop.r2 = 0.025; % outer electrode radius [m]
prop.r1 = 0.024; % inner electrode radius [m]
prop.L = 0.1;    % length of APM [m]
RPM = 13350; % rotational speed [rpm]
prop.omega = RPM*2*pi/60; % rotational speed [rad/s]
prop.Q = 1.02E-3/60; % aerosol flowrate [m^3/s]
prop.T = 298; % system temperature [K]
prop.p = 1; % system pressure [atm]

end


