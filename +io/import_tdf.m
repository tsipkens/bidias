
% IMPORT_TDF  Import general tandem data files. 
%  
%  Cambustion file format containing 
%  
%  AUTHOR: Timothy Sipkens, 2023-09-22

function [data, dm, sp, prop1, prop2] = import_tdf(fname)

addpath tfer_pma;  % in case not already added

opts = detectImportOptions(fname, 'FileType', 'text');
read = readmatrix(fname, 'FileType', 'text');

in = readtable(fname);

head = readcell(fname);
head = head(1:opts.VariableNamesLine-1, :);

class1 = head(2:3, :);
class2 = head(4:5, 1:end-1);


% Assume CPMA, then DMA...

%== CPMA read =============================================%
idx = find(contains(class1(1,:), 'Sample flow'));
Qsmpl = class1{2, idx};

prop1 = kernel.prop_pma;
prop1 = prop_update_flow(prop1, Qsmpl/1000/60);
prop1.T = mean(in.Temperature_C_) + 273;  %  use average for now
prop1.p = mean(in.Pressure_Pa_) ./ 101325;

sp = get_setpoint(prop1,...
    'omega', in.Speed_rad_s_, 'V', in.Voltage_V_); % get CPMA setpoint information


%== DMA read ==============================================%
idx = find(contains(class2(1,:), 'Sample flow'));
Qsmpl = class2{2, idx};

idx = find(contains(class2(1,:), 'Sheath flow'));
Qsheath = class2{2, idx};

opts = struct();
opts.params = 'custom';
opts.prop.Q_s = Qsmpl/60/1000;  % sample flow [m^3/s]
opts.prop.Q_a = Qsmpl/60/1000;  % aerosol flow [m^3/s]
opts.prop.Q_c = Qsheath/60/1000;  % sheath flow [m^3/s]
opts.prop.Q_m = Qsheath/60/1000;  % exhaust flow [m^3/s]
opts.prop.T = mean(in.Temperature_C_2);
opts.prop.p = mean(in.Pressure_kPa_2) ./ 101.325;
opts.prop.L = 0.44369;  % length of chamber [m]
opts.prop.R2 = 0.00937; % outer electrode radius [m]
opts.prop.R1 = 0.01961; % inner electrode radius [m]
prop2 = kernel.prop_dma(opts);  % fill out rest of the parameters

dm = in.Dm_nm_2;

data = in.Conc;

end
