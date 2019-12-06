
% GRID2SP Convert a grid used for data to a structured array of PMA setpoints.
% Author: Timothy Sipkens, 2019-12-06
%=========================================================================%

function [sp] = grid2sp(prop,grid_b,varargin)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = []; end
if isempty(prop); prop = kernel.prop_pma; end
    % import properties of PMA
    % use default properties selected by prop_pma function

if isempty(varargin); varargin = {'Rm',3}; end
    % if second component of setpoint is not specified
    % use a resolution of Rm = 3
%-------------------------------------------------------------------------%


m_star = grid_b.elements(:,1);

for ii=length(m_star):-1:1
    sp(ii) = tfer_pma.get_setpoint(prop,...
        'm_star',m_star(ii),varargin{:});
end

end

