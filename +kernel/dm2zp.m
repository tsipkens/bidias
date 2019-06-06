
function [B,Zp] = dm2zp(d,z,T,p)
% DP2ZP Calculate electric mobility from a vector of mobility diameter.
% Author: Timothy Sipkens, 2019-1-2
%
% Note: Combines functions from Buckley et al. and Olfert laboratory.


%-- Parse inputs ----------------------------------------%
if ~exist('z','var')
    z = 1; % include diffusion by default
elseif isempty(z)
    z = 1;
end


%-- Perform calculation ---------------------------------%
e = 1.6022e-19; % electron charge [C]
if nargin<=2 % if pressure and temperature are not specified (Buckley et al. / Davies)
    mu = 1.82e-5; % gas viscosity [Pa*s]
    B = Cc(d)./(3*pi*mu.*d); % mechanical mobility
    
else % from Olfert laboratory / Kim et al.
    S = 110.4; % temperature [K]
    T_0 = 296.15; % reference temperature [K]
    vis_23 = 1.83245*10^-5; % reference viscosity [kg/(m*s)] or [Pa*s] or [N*s/m2]
    mu = vis_23*((T/T_0)^1.5)*((T_0+S)/(T+S)); % gas viscosity
        % Kim et al. (2005), ISO 15900 Eqn 3
        
    B = Cc(d,p,T)./(3*pi*mu.*d); % mechanical mobility
    
end

Zp = B.*e.*z;

end


function out = Cc(d,p,T)
% CC Cunningham slip correction factor.

if nargin==1 % if pressure and temperature are not specified (Buckley et al. / Davies)
    mfp = 66.5e-9; % mean free path
    
    % for air, from Davies, 1945
    A1 = 1.257;
    A2 = 0.4;
    A3 = 0.55;
    
else % from Olfert laboratory / Kim et al.
    S = 110.4; % temperature in K
    mfp_0 = 6.730*10^-8; % mean free path of gas molecules in air at reference conditions [m]
    T_0 = 296.15; % reference temperature [K]
    p_0 = 101325; % reference pressure, [Pa] (760 mmHg to Pa)
    
    p = p*p_0;
    
    % Kim et al. (2005) (10.6028/jres.110.005), ISO 15900 Eqn 4
    mfp = mfp_0*(T/T_0)^2*(p_0/p)*((T_0+S)/(T+S)); % mean free path
    
    A1 = 1.165;
    A2 = 0.483;
    A3 = 0.997/2;
    
end

Kn = (2*mfp)./d; % Knudsen number
out = 1 + Kn.*(A1 + A2.*exp(-(2*A3)./Kn)); % Cunningham slip correction factor

end
