
% IMPORT_TDF  Import general tandem data files. 
%  
%  TDF are tab-delimited files that contain particle counts alongside a
%  large range of classifier information. 
%  
%  AUTHOR: Timothy Sipkens, 2022-12-13

function [data, dm, sp, prop1, prop2] = import_tdf(fname)

addpath tfer_pma;  % in case not already added

opts = detectImportOptions(fname, 'FileType', 'text');
read = readmatrix(fname, 'FileType', 'text');

% Read components of header.
fid = fopen(fname);
tt = 1;
tline = fgetl(fid);
while ischar(tline)
    if tt == 4
        hline = split(tline, '	');
        l1 = str2num(hline{3});  % read number of reported data columns
        Q1 = str2num(hline{5});
    end

    if tt == 8
        hline = split(tline, '	');
        l2 = str2num(hline{3});  % read number of reported data columns
        Q2 = str2num(hline{4});
    end
    
    tt = tt + 1;
    tline = fgetl(fid);
end
fclose(fid);

% For CPMA.
s1 = read(:,1);
dim1 = length(s1);
w = read(:, 2);
V = read(:, 3);

% Setup prop for the PMA.
prop1 = prop_pma('tfer_pma/prop/olfert');
prop1 = prop_update_flow(prop1, Q1 / 1000 / 60);
prop1 = prop_update_massmob(prop1, 'Dm', 2.48, 'rho100', 510);
sp = get_setpoint(prop1, 'V', V, 'omega', w);

%-- DMA properties ---------%
dim2 = (size(read, 2) - l1) / l2;
dm = read(:, l1 + 1:l2:end);
dm = dm(1,:);

Qsh = mean(mean(read(:, l1 + 2:l2:end))) / 1000 / 60;
prop2 = kernel.prop_dma();
prop2.Q_c = Qsh;
prop2.Q_m = Qsh;
prop2.Q_s = Q2 / 1000 / 60;
prop2.Q_a = Q2 / 1000 / 60;

data = read(:, l1 + l2:l2:end);

end
