
% GEN_GD  Generate
% Author: Timothy Sipkens, 2018-10-22
% 
% Two possible options for  will be a series of name value pairs, including:
%   'l1' or 's1' - Correlation length for variable 1 (e.g. mass).
%   'l2' or 's2' - Correlation length for variable 2 (e.g. mobility diameter).
%   'R12' - Correlation between variables.
%   'Dm' - Mass-mobility exponent or slope between variables
%=========================================================================%

function [Gd, sg] = gen_Gd(varargin)

sg = struct('l1',[],'l2',[],...
    'R12',[],'Dm',[]); % structure of Gd properties


if length(varargin)==1 % Gd was supplied, return correlation lengths, corr., etc.
    Gd = varargin{1};
    
    
else % otherwise, name-value pairs supplied
    for ii=1:2:length(varargin) % loop through name-values pairs
        % allow for 's1' or 's2' in the place of 'l1' and 'l2'
        if strcmp(varargin{ii},'s1'); varargin = 'l1'; end
        if strcmp(varargin{ii},'s2'); varargin = 'l2'; end

        sg.(varargin{ii}) = varargin{ii+1}; % copy Gd properties to structure
    end
    
    
    Gd = [sg.l1^2, 0; ...
        0, sg.l2^2]; % initiate Gd using correlation lengths
    
    
    %-- Incorporate correlation information ------------------------------%
    if ~isempty(sg.R12) % OPTION 1: By default, use supplied value of R12
        Gd(1,2) = sg.R12 * sg.l1 * sg.l2;
    
    else % OPTION 2: Alternatively, use mass-mobility exponent / slope
        Gd(1,2) = sg.Dm * sg.l2;
    
    end
    Gd(2,1) = Gd(1,2);
    %---------------------------------------------------------------------%
end


%-- Update sg structure --------------------------------------------------%
sg.l1 = sqrt(Gd(1,1));
sg.l2 = sqrt(Gd(2,2));
sg.R12 = Gd(1,2) / (sg.l1 * sg.l2);
sg.Dm = Gd(1,2)/Gd(2,2);
%-------------------------------------------------------------------------%


end

