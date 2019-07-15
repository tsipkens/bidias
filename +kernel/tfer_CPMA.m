
function [Lambda,prop] = tfer_CPMA(m_star,m,d,z,prop,opt,varargin)
% TFER_CPMA     Bridging function used to evaluate CPMA transfer function.
% Author:       Timothy Sipkens, 2018-12-27
%
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
%
% Note:
%   Some of this code is adapted from Buckley et al. (2017) and Olfert
%   laboratory.
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = []; end
if ~exist('opt','var'); opt = []; end
if ~exist('varargin','var'); varargin = []; end

if isempty(prop); prop = kernel.prop_CPMA; end % import properties of CPMA
if isempty(opt); opt = 'B_diff'; end
    % by default use finite difference solution
if isempty(varargin); varargin = {'Rm',3}; end
    % by default use resolution of 3


addpath('UBC-tfer-PMA'); % include PMA transfer functions as submodule,
                         % included as the UBC-tfer-PMA directory

fun = str2func(['tfer_PMA.tfer_',opt]); % call relevant function from submodule
Lambda = fun(m_star,m,d,z,prop,varargin{1},varargin{2})'; % CPMA transfer function


end
