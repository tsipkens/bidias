
% IMPORT_C Read tandem CPMA-SMPS data from a list-format (newer format) CSV file.
% Author: Timothy Sipkens, 2019-12-17
% 
% Note:
%   The prop_pma input would override deriving some of the properties from
%   the imported data file.
%=========================================================================%

function [data0,d_star0,sp0,prop_dma,prop_pma,out] = ...
    import_c(fnames,prop_pma)

%-- Parse inputs ---------------------------------------------------------%
addpath tfer_pma; % add tfer_pma dependency

% check if input is a filename structure, or convert to structure
if ~isstruct(fnames)
    t0 = fnames;
    fnames = struct();
    [fnames.folder,fnames.name,ext] = fileparts(t0);
    fnames.name = [fnames.name,ext];
end

% if not supplied get properties structure
if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma(' CPMA'); end
%-------------------------------------------------------------------------%


% initiate diameters
data0 = [];
d_star0 = [];


disp('Reading files...');
N = length(fnames); % number of files to read
if N>1; tools.textbar(0); end % initiate commandline textbar
for ff=1:N % loop through files
    fn = [fnames(ff).folder,'\',fnames(ff).name];
        % ffth file name
    
    
    opts = detectImportOptions(fn,'NumHeaderLines',1);
        % default options
    
    %== Proceed with reading data ================%
    tab0 = readtable(fn, opts);
    
    %-- PMA setpoints -----------------------------------------------%
    m_star = tab0.m_fg_;
    V = tab0.V_V__1;
    omega = tab0.w_rad_s_; % centerline radial speed
    Rm = []; % CPMA implied resolution, currently N/A
    p = tab0.P_Pa_; % PMA pressure
    T = tab0.T_C_; % PMA tempreature
    
    
    
    prop_pma.Q = mean(tab0.Qa_LPM__1) / 1000 / 60; % current assumed constant over all setpoint
            % and does not support changing Q during measurements
    prop_pma.v_bar = prop_pma.Q/prop_pma.A; % average flow velocity
    
    clear sp; % reset PMA setpoints
    prop_pma.T = mean(T) + 273.15; % pressure converted to Kelvin
    prop_pma.p = mean(p) / 101325; % pressure converted to atm
    % sp(ii,1) = tfer_pma.get_setpoint(prop_pma,...
    %     'm_star',m_star(ii).*1e-18,'Rm',Rm(ii));
    sp = get_setpoint(prop_pma,...
        'V', V, 'omega', omega);
    %----------------------------------------------------------------%


    %-- DMA setpoints -----------------------------------------------%
    d_star = tab0.dm_nm_; % DMA setpoints
    prop_dma = kernel.prop_dma(' Electrostatic Classifier Model 3080');
    prop_dma.Q_a = prop_pma.Q; % aerosol flow (same as PMA)
    prop_dma.Q_s = prop_dma.Q_a;
    prop_dma.Q_c = mean(tab0.Qsh_LPM_) / 1000 / 60; % shield flow
    prop_dma.Q_m = prop_dma.Q_c;
    %----------------------------------------------------------------%
    
    
    data = tab0.N_p_cc_; % data as a vector
    
    
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
