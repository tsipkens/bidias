
% GEN_GD  Generate the Gd matrix from a series of name-value pairs.
% Author: Timothy Sipkens, 2020-07-17
% 
% Inputs: 
%   OPTION 1: One can enter a series of name-value pairs.
%   At least, three of the following arguments must be supplied:
%       'l1' or 's1' - Correlation length for variable 1 (e.g. mass).
%       'l2' or 's2' - Correlation length for variable 2 (e.g. mobility diameter).
%       'R12' - Correlation between variables.
%       'Dm' - Mass-mobility exponent or slope between variables
%         e.g., varargin = {'l1',0.2,'Dm',3,'R12',0.99}
%   
%   OPTION 2: Alternatively, one can supply Gd, which is used to fill out the 
%   sg structure, containing the above properties. Remaining name-value
%   pairs are used to output a new Gd by modifying l2.
%         e.g., varargin = {Gd,'R12',0.99}
%=========================================================================%

function [Gd, sg] = gen_Gd(varargin)

sg = struct('l1',[],'l2',[],...
    'R12',[],'Dm',[]); % structure of Gd properties


%== OPTION 2: Gd was supplied, return correlation lengths, corr., etc. ===%
%   Update parameters if additional arguments are supplied.
if isa(varargin{1}, 'double')
    Gd = varargin{1};
    
    % evaluate parameter set implied by Gd, before modifiction
    sg.l1 = sqrt(Gd(1,1));
    sg.l2 = sqrt(Gd(2,2));
    sg.R12 = Gd(1,2) / (sg.l1 * sg.l2);
    sg.Dm = Gd(1,2) / Gd(2,2);
    
    % loop through name-values pairs
    for ii=2:2:length(varargin)
        if strcmp(varargin{ii},'s1'); varargin{ii} = 'l1'; end % allow for 's1' in the place of 'l1'
        if or(strcmp(varargin{ii},'s2'),strcmp(varargin{ii},'12')) % cannot use l2
            warning
            continue;
        end
        sg.(varargin{ii}) = varargin{ii+1}; % copy Gd properties to structure
    end
    sg.l2 = sg.R12 / sg.Dm * sg.l1;
    

%== OPTION 1: Otherwise, name-value pairs supplied =======================%
else
    % loop through name-values pairs
    for ii=1:2:length(varargin) 
        % allow for 's1' or 's2' in the place of 'l1' and 'l2'
        if strcmp(varargin{ii},'s1'); varargin{ii} = 'l1'; end
        if strcmp(varargin{ii},'s2'); varargin{ii} = 'l2'; end
        sg.(varargin{ii}) = varargin{ii+1}; % copy Gd properties to structure
    end
    
    % update correlation lengths, if necessary
    if isempty(sg.l1) % use Dm, R12, and l2 if l1 is missing
        sg.l1 = sg.Dm / sg.R12 * sg.l2;
    elseif isempty(sg.l2) % use Dm, R12, and l1 if l2 is missing
        sg.l2 = sg.R12 / sg.Dm * sg.l1;
    end
end


% initiate Gd using correlation lengths
% re-evaluates to Gd if only Gd was supplied as input to method
Gd = [sg.l1^2, 0; ...
    0, sg.l2^2];


% incorporate correlation lengths
if ~isempty(sg.R12) % OPTION 1: By default, use supplied value of R12
    Gd(1,2) = sg.R12 * sg.l1 * sg.l2;

else % OPTION 1: Alternatively, use mass-mobility exponent / slope
    Gd(1,2) = sg.Dm * sg.l2 .^ 2;

end
Gd(2,1) = Gd(1,2);


%-- Update sg structure --------------------------------------------------%
sg.l1 = sqrt(Gd(1,1));
sg.l2 = sqrt(Gd(2,2));
sg.R12 = Gd(1,2) / (sg.l1 * sg.l2);
sg.Dm = Gd(1,2) / Gd(2,2);
%-------------------------------------------------------------------------%


% return error if correlation is unphysical
if sg.R12>=1
    error('Correlation implied by Gd matrix exceeds unity.');
end

end

