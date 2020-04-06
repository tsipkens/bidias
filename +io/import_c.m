
% IMPORT_C Read tandem CPMA-SMPS data from a list-format (newer format) CSV file.
% Author: Timothy Sipkens, 2019-12-17
% 
% Note:
%   The prop_pma input would override deriving some of the properties from
%   the imported data file.
%=========================================================================%

function [data0,d_star0,sp0,prop_dma,prop_pma,out] = ...
    import_c(fnames,prop_pma)

if ~isstruct(fnames)
    t0 = fnames;
    fnames = struct();
    [fnames.folder,fnames.name,ext] = fileparts(t0);
    fnames.name = [fnames.name,ext];
end

data0 = [];
d_star0 = [];

disp('Reading files...');

N = length(fnames); % number of files to read
if N>1; tools.textbar(0); end
for ff=1:N
    fn = [fnames(ff).folder,'\',fnames(ff).name];
    
    
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

    if ~exist('prop_pma','var'); prop_pma = []; end
    if isempty(prop_pma); prop_pma = kernel.prop_pma(' CPMA'); end
    
    prop_pma.Q = mean(table2array(t(:,9)))/1000/60; % current averages over all setpoint
            % and does not support changing Q during measurements
    prop_pma.v_bar = prop_pma.Q/prop_pma.A; % average flow velocity
    
    clear sp; % reset PMA setpoints
    for ii=1:length(V)
        prop_pma.T = T(ii)+273.15; % pressure converted to Kelvin
        prop_pma.p = p(ii)/101325; % pressure converted to atm
        % sp(ii,1) = tfer_pma.get_setpoint(prop_pma,...
        %     'm_star',m_star(ii).*1e-18,'Rm',Rm(ii));
        sp(ii,1) = tfer_pma.get_setpoint(prop_pma,...
            'V',V(ii),'omega',omega(ii));
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
