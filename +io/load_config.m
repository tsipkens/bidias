
% LOAD_CONFIG  Loads settings from configuration file (YML). 
%  Files are loaded in order supplied, overwriting properties where
%  relevant. 
%  
%  AUTHOR: Timothy Sipkens, 2021-03-25

function config = load_config(fnames)

if ~iscell(fnames); fnames = {fnames}; end

config = struct();
for ii=1:length(fnames)
    
    config0 = io.read_yml(fnames{ii});  % read new settings
    
    % Copy (or overwrite) existing settings.
    f = fieldnames(config0);
    for jj = 1:length(f)
        config.(f{jj}) = config0.(f{jj});
    end
    
end

end
