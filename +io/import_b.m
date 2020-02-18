
% IMPORT_B Read tandem CPMA-SMPS data from a list-format (Tyler format) CSV file.
% Author: Timothy Sipkens, 2019-12-17
%=========================================================================%

function [data0,d_star0,sp0,prop_dma,prop_pma,out] = ...
    import_b(fns)

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
    [~] = split(fgets(fid),',');
    line2 = split(fgets(fid),',');
    line3 = split(fgets(fid),',');
    
    ii = 3; % automatically detect size of header
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
    p = table2array(t(:,6)); % PMA pressure
    T = table2array(t(:,7)); % PMA tempreature

    prop_pma = kernel.prop_pma(line2{2});
    prop_pma.mass_mob_pref = 0.0612; % assume CPMA uses soot properties
    prop_pma.mass_mob_exp = 2.48;
    
    clear sp; % reset PMA setpoints
    for ii=1:length(V)
        prop_pma.T = T(ii)+273.15; % pressure converted to Kelvin
        prop_pma.p = p(ii)/101325; % pressure converted to atm
        sp(ii,1) = tfer_pma.get_setpoint(prop_pma,...
            'omega',omega(ii),'V',V(ii));
    end
    %----------------------------------------------------------------%


    %-- DMA setpoints -----------------------------------------------%
    d_star = table2array(t(:,11)); % DMA setpoints
    prop_dma = kernel.prop_dma(line3{3});
    %----------------------------------------------------------------%


    data = table2array(t(:,19)); % data as a vector
    
    
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