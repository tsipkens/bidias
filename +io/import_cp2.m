
% IMPORT_CP2  Importer for CPMA-SP2 data. 
%  
%  AUTHOR: Timothy Sipkens, 2023-10-23

function [R, sp, prop, in, t_span0, dt] = import_cp2(fn1, fn2, shift)

if ~exist('shift', 'var'); shift = []; end
if isempty(shift)
    shift = seconds(0);  % default is no time mistmatch
end

tools.textheader('Importing tandem CPMA-SP2 data');

% Import CPOF file.
disp(' Reading CPOF file...');
[sp, prop, ~, t_span0] = io.read_cpof(fn1);
t_span0 = t_span0 + shift;  % shift datetimes, if necessary
tools.textdone();

% Read in SP2 data file.
disp(' Reading SP2 data file...');
in = readtable(fn2);
tools.textdone();

% Extract time span from CPOF file. 
% Assume only single time at each setpoint.
% Currently also do not account for buffer time between scans. 
R = {};  % response
t_span(length(sp), 2) = datetime();
for ii=1:length(sp)
    
    % Find relevant SP2 data. 
    fl_sp2 = and(in.TimeDate > t_span0(ii,1), ...
                 in.TimeDate < t_span0(ii,2));
    
    % Get actual time span in SP2 data.
    if sum(fl_sp2) ~= 0
        t_span(ii,1) = min(in.TimeDate(fl_sp2));
        t_span(ii,2) = max(in.TimeDate(fl_sp2));
    end
    
    R{ii} = in.BHBL_BCmass(fl_sp2);

end

dt = t_span(:,2) - t_span(:,1);

tools.textheader();

end
