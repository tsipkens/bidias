
function [B,Zp] = dm2zp(d,z,T,p)
% DP2ZP Calculate electric mobility from a vector of mobility diameter.
% Adapted from: Buckley et al. and Olfert laboratory.
% Author: Timothy Sipkens, 2019-01-02
% Note: Temperature and pressure must both be specified before they are
%   used.
%
%-------------------------------------------------------------------------%
% Inputs:
%   d           Particle mobility diameter
%   z           Integer charge state
%   T           System temperature  
%                   (Optional, if not specified use code from Buckley et al.)
%   P           System pressure     (Optional, same as above)
%
% Outputs:
%   Zp          Electromobility
%   B           Mechanical mobility
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('z','var') % if integer charge state not specified use z = 1
    z = 1;
elseif isempty(z)
    z = 1;
end


%-- Perform calculation --------------------------------------------------%
e = 1.6022e-19; % electron charge [C]
if nargin<=3 % if pressure and temperature are not specified (Buckley et al. / Davies)
    mu = 1.82e-5; % gas viscosity [Pa*s]
    B = Cc(d)./(3*pi*mu.*d); % mechanical mobility
    
else % from Olfert laboratory / Kim et al.
    S = 110.4; % temperature [K]
    T_0 = 296.15; % reference temperature [K]
    vis_23 = 1.83245*10^-5; % reference viscosity [kg/(m*s)] or [Pa*s] or [N*s/m2]
    mu = vis_23*((T/T_0)^1.5)*((T_0+S)/(T+S)); % gas viscosity
        % Kim et al. (2005), ISO 15900 Eqn 3
        
    B = Cc(d,T,p)./(3*pi*mu.*d); % mechanical mobility
    
end

Zp = B.*e.*z;

end


function Cc = Cc(d,T,p)
% CC Function to evaluate Cunningham slip correction factor.
% Author: Timothy Sipkens, 2019-01-02
%
%-------------------------------------------------------------------------%
% Inputs:
%   d           Particle mobility diameter
%   T           System temperature  
%                   (Optional, if not specified use code from Buckley et al.)
%   P           System pressure     (Optional, same as above)
%
% Outputs:
%   Zp          Electromobility
%   B           Mechanical mobility
%-------------------------------------------------------------------------%

if nargin==1 % if pressure and temperature are not specified (Buckley et al. / Davies)
    mfp = 66.5e-9; % mean free path
    
    % for air, from Davies (1945)
    A1 = 1.257;
    A2 = 0.4;
    A3 = 0.55;
    
else % from Olfert laboratory / Kim et al.
    S = 110.4; % temperature in K
    mfp_0 = 6.730*10^-8; % mean free path of gas molecules in air at reference conditions [m]
    T_0 = 296.15; % reference temperature [K]
    p_0 = 101325; % reference pressure, [Pa] (760 mmHg to Pa)
    
    p = p*p_0;
    
    % Kim et al. (2005) (doi:10.6028/jres.110.005), ISO 15900 Eqn 4
    mfp = mfp_0*(T/T_0)^2*(p_0/p)*((T_0+S)/(T+S)); % mean free path
    
    A1 = 1.165;
    A2 = 0.483;
    A3 = 0.997/2;
    
end

Kn = (2*mfp)./d; % Knudsen number
Cc = 1 + Kn.*(A1 + A2.*exp(-(2*A3)./Kn)); % Cunningham slip correction factor

end
