
% TFER_PMA     Bridging function used to evaluate particle mass analyer (PMA) transfer function.
% Author:       Timothy Sipkens, 2018-12-27
%-------------------------------------------------------------------------%
% Inputs:
%   sp          Structure defining various setpoint parameters 
%               (e.g. m_star, V). Use 'get_setpoint' method to generate 
%               this structure.
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length) (optional)
%   opt         Alphanumeric code for transfer function evaluation method
%               (optional)
%
% Outputs:
%   Lambda      Transfer function
%   prop        CPMA device settings
%=========================================================================%

function [Lambda,prop] = tfer_pma(sp,m,d,z,prop,opt)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('opt','var'); opt = []; end

if isempty(opt); opt = '1C_diff'; end
    % by default, use Taylor series solution baout rc (Case 1C) with diffusion
    % see Sipkens et al., Aerosol Sci. Technol. (2019) for more information
%-------------------------------------------------------------------------%


fun = str2func(['tfer_pma.tfer_',opt]); % call relevant function from submodule
Lambda = fun(sp,m,d,z,prop)'; % CPMA transfer function


end
