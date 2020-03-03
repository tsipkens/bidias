
% IMPORT_C Read tandem CPMA-SMPS data from a list-format (newer format) CSV file.
% Author: Timothy Sipkens, 2019-12-17
%=========================================================================%

function [data0,d_star0,sp0,prop_dma,prop_pma,out] = ...
    import_c(fns)

if ~isstruct(fns)
    t0 = fns;
    fns = struct();
    [fns.folder,fns.name,ext] = fileparts(t0);
    fns.name = [fns.name,ext];
end

data0 = [];
d_star0 = [];

disp('Reading files...');

N = length(fns); % number of files to read
if N>1; tools.textbar(0); end
for ff=1:N
    fn = [fns(ff).folder,'\',fns(ff).name];
    
    
    %== Open file and read data =========%
    fid = fopen(fn);
    
    ii = 0; % automatically detect size of header
    line = {''};
    while ~strcmp(line{1},'Point')
        line = split(fgets(fid),',');
        ii = ii+1;
    end
    
    [~] = fclose(fid);
    
    
    %== Proceed with reading data ========%
    t = readtable(fn, 'HeaderLines',ii-1);
    
    %-- PMA setpoints -----------------------------------------------%
    m_star = table2array(t(:,3));
    V = table2array(t(:,4));
    omega = table2array(t(:,5)); % centerline radial speed
    Rm = table2array(t(:,6)); % CPMA implied resolution
    p = table2array(t(:,7)); % PMA pressure
    T = table2array(t(:,8)); % PMA tempreature

    prop_pma = kernel.prop_pma(' CPMA');
    prop_pma.mass_mob_pref = 0.0612; % assume CPMA uses soot properties
    prop_pma.mass_mob_exp = 2.48;
    % prop.mass_mob_pref = 524;
    % prop.mass_mob_exp = 3;
    
    prop_pma.Q = mean(table2array(t(:,9)))/1000/60; % current averages over all setpoint
            % and does not support changing Q during measurements
    
    clear sp; % reset PMA setpoints
    for ii=1:length(V)
        prop_pma.T = T(ii)+273.15; % pressure converted to Kelvin
        prop_pma.p = p(ii)/101325; % pressure converted to atm
        sp(ii,1) = tfer_pma.get_setpoint(prop_pma,...
            'm_star',m_star(ii).*1e-18,'Rm',Rm(ii));
    end
    %----------------------------------------------------------------%


    %-- DMA setpoints -----------------------------------------------%
    d_star = table2array(t(:,12)); % DMA setpoints
    prop_dma = kernel.prop_dma(' Electrostatic Classifier Model 3080');
    prop_dma.Q_a = mean(table2array(t(:,9)))/1000/60; % aerosol flow (same as PMA)
    prop_dma.Q_s = prop_dma.Q_a;
    prop_dma.Q_c = mean(table2array(t(:,14)))/1000/60; % shield flow
    prop_dma.Q_m = prop_dma.Q_c;
    %----------------------------------------------------------------%
    
    
    data = table2array(t(:,20)); % data as a vector
    
    
    %-- Append to arrays --------------------------------------------%
    if ~exist('sp0','var')
        sp0 = sp;
    else
        sp0 = [sp0;sp];
    end
    data0 = [data0;data];
    d_star0 = [d_star0;d_star];
    out{ff} = [d_star,m_star,data];
    %----------------------------------------------------------------%
    
    
    if N>1; tools.textbar(ff./N); end
end

disp('Complete.');
disp(' ');


end
