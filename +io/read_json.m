
% READ_JSON  Read JSON structured configuration files. 
%  Allows for C++ or Javscript style commenting.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-20

function results = read_json(file)

fid = fopen(file);
raw = fread(fid, inf);  % raw file contents
str = char(raw');  % transform to char
fclose(fid);

% Remove comments.
str = erase(erase(eraseBetween( ...
    erase(eraseBetween(str, "//", newline), "//"), ...
    "/*", "*/"), "/*"), "*/");

results = jsondecode(str);

% Attempt to interpret Matlab expressions.
f_results = fields(results);
for ii=1:length(f_results)
    t0 = results.(f_results{ii});
    
    if isa(t0, 'char')
        [converted, success] = str2num(t0);
        if success
            results.(f_results{ii}) = converted;
        end
    end
end

end

