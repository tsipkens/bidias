
% TFER_PMA     Bridging function used to evaluate particle mass analyer (PMA) transfer function.
% Author:       Timothy Sipkens, 2018-12-27
%=========================================================================%

function [Lambda,prop] = tfer_pma(m_star,m,d,z,prop,opt,varargin)
%--------------------------------------------------------------------------%
% Inputs:
%   m_star      Mass corresponding to the measurement set point of the APM
%   d           Particle mobility diameter, can be vector [nm]
%   m           Particle mass, can be vector of same length as d
%   z           Integer charge state, scalar
%   prop        CPMA device settings        (Optional)
%   opts        Case used for tfer. func. evaluation
%                   (Optional, default is Case B w/ diffusion)
%   varargin    Name-value pairs for setpoint
%                   (Optional, default {Rm, 3})
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Outputs:
%   Lambda      CPMA transfer function
%   prop        CPMA device settings
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = []; end
if ~exist('opt','var'); opt = []; end
if ~exist('varargin','var'); varargin = []; end

if isempty(prop); prop = kernel.prop_pma; end
    % import properties of PMA
    % use default properties selected by prop_pma function

if isempty(opt); opt = '1C_diff'; end
    % by default, use Taylor series solution baout rc (Case 1C) with diffusion
    % see Sipkens et al., Aerosol Sci. Technol. (2019) for more information

if isempty(varargin); varargin = {'Rm',3}; end
    % by default, use resolution of 3
%-------------------------------------------------------------------------%


fun = str2func(['tfer_pma.tfer_',opt]); % call relevant function from submodule
Lambda = fun(m_star,m,d,z,prop,varargin{1},varargin{2})'; % CPMA transfer function


end
