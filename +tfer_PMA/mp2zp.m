function [Zp,B,d] = mp2zp(m,z,T,P)
% MP2ZP Calculate electric mobility from a vector of particle mass.
% Author: Timothy Sipkens, 2019-1-2
% Note: Uses mass-mobility relationship to first convert to a mobility
%   diameter and then estimates the mobility using dm2zp.
% 
%-------------------------------------------------------------------------%
% Inputs:
%   m           Particle mass
%   z           Integer charge state
%   T           System temperature
%   P           System pressure
%
% Outputs:
%   Zp          Electromobility
%   B           Mechanical mobility
%   d           Mobility diameter (implied by mass-mobility relation)
%-------------------------------------------------------------------------%


%-- Invoke mass-mobility relation ----------------------------------------%
mass_mob_pref = 524;
mass_mob_exp = 3;
d = (m./mass_mob_pref).^(1/mass_mob_exp);
    % use mass-mobility relationship to get mobility diameter

    
%-- Use mobility diameter to get particle electro and mechanical mobl. ---%
if nargin<3
    [Zp,B] = tfer_PMA.dm2zp(d,z);
else
    [Zp,B] = tfer_PMA.dm2zp(d,z,T,P);
end

end

