
% LOAD_CONFIG  Loads settings from configuration file (YML, YAML, or JSON). 
%  Files are loaded in order supplied, overwriting properties where
%  relevant. 
%  
%  AUTHOR: Timothy Sipkens, 2021-03-25

function config = load_config(fnames)

if ~iscell(fnames); fnames = {fnames}; end

config = struct();
for ii=1:length(fnames)
    
    if strcmp(fnames{ii}(end-3:end), 'json')
        config0 = io.read_json(fnames{ii});  % read new settings
        
    elseif or(strcmp(fnames{ii}(end-2:end), 'yml'), strcmp(fnames{ii}(end-3:end), 'yaml'))
        config0 = io.read_yml(fnames{ii});  % read new settings
        
    else; error('Failed to load config file.'); end
    
    % Copy (or overwrite) existing settings.
    f = fieldnames(config0);
    for jj = 1:length(f)
        config.(f{jj}) = config0.(f{jj});
    end
    
end

end
