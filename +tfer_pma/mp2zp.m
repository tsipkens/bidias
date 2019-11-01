
% MP2ZP     Calculate electric mobility from a vector of particle mass.
% Author:   Timothy Sipkens, 2019-01-02
%=========================================================================%

function [Zp,B,d] = mp2zp(m,z,T,P,prop)

%-------------------------------------------------------------------------%
% Inputs:
%   m           Particle mass
%   z           Integer charge state
%   T           System temperature
%   P           System pressure
%   prop        CPMA/DMA properties structure
%
% Outputs:
%   Zp          Electromobility
%   B           Mechanical mobility
%   d           Mobility diameter (implied by mass-mobility relation)
% 
% Note:
%   Uses mass-mobility relationship to first convert to a mobility
%   diameter and then estimates the mobility using dm2zp.
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('T','var'); T = []; end
if ~exist('P','var'); P = []; end

if ~exist('prop','var'); prop = []; end
if or(isempty(prop),...
        ~and(isfield(prop,'mass_mob_pref'),...
        isfield(prop,'mass_mob_exp'))) % get parameters for the mass-mobility relation
    mass_mob_pref = 0.0612;
    mass_mob_exp = 2.48;
else
    mass_mob_pref = prop.mass_mob_pref;
    mass_mob_exp = prop.mass_mob_exp;
end
% if isempty(prop); mass_mob_pref = 524; mass_mob_exp = 3; end
%-------------------------------------------------------------------------%


d = (m./mass_mob_pref).^(1/mass_mob_exp);
    % use mass-mobility relationship to get mobility diameter

    
%-- Use mobility diameter to get particle electro and mechanical mobl. ---%
if or(isempty(T),isempty(P))
    [Zp,B] = tfer_pma.dm2zp(d,z);
else
    [Zp,B] = tfer_pma.dm2zp(d,z,T,P);
end

end

