
% LINECOUNT   Count lines in the file.
% Author: Timothy Sipkens, 2020-04-12
%=========================================================================%

function n = linecount(fname)

fid = fopen(fname);

n = 0;
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;
end

fclose(fid);

end