
% IMPORT_2  Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, 2019-12-17
%
% Note: fn_cpma can be either a single file name or a cell array of file
% names. The times are combined / stacked vertically before a search for
% the times corresponding to the DMA time stamps is performed. CPMA data
% files should be given in chronological order.
%=========================================================================%

function [data,d_star,sp,prop_dma,prop_pma] = import_2(fn_smps,fn_cpma,prop_pma)


%== Read SMPS/CPC file =========================================%
[data,d_star,prop_dma,time_dma] = io.import_dma(fn_smps);
%===============================================================%



%== Read CPMA file =============================================%
addpath('tfer_pma');

if ~iscell(fn_cpma); fn_cpma = {fn_cpma}; end

opts = detectImportOptions(fn_cpma{1}); % default options
for ff=1:length(fn_cpma)
    if ~exist('time_cpma','var'); idx0 = 1;
    else; idx0 = [idx0, length(time_cpma)]; end
    
    t0 = readtable(fn_cpma{ff}, opts); % read CPMA file
    time_cpma0 = datetime(t0.Date_Time); % get times to compare against DMA data
    
    if ~exist('t1','var'); t1 = t0; time_cpma = time_cpma0;
    else; t1 = [t1;t0]; time_cpma = [time_cpma;time_cpma0]; end
end


[gr1,gr2] = ndgrid(time_dma, time_cpma); % grid of CPMA/SMPS times
[~,idx_t] = max(gr1<gr2,[],2); % get first CPMA time corresponding to DMA scan

if any(idx_t==idx0)
    warning(['CPMA data files do not contain correct times. ' ...
        'Only DMA data output.']);
    sp = get_setpoint();
    return;
end


% voltage and speed used to define CPMA setpoint
omega = t1.ClassSpFB_rad_s__(idx_t); % speed of CPMA (no averaging)
V = t1.AdjVoltageFB_V__(idx_t); % voltage of CPMA (no averaging)


% alternate voltages / speeds (unused)
V_raw = t1.RawVoltageFB_V__(idx_t);
omega1 = t1.InnerSpFB_rad_s__(idx_t);
omega2 = t1.OuterSpFB_rad_s__(idx_t);


% if prop_pma is not specified
if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end


prop_pma.T = mean(t1.Temperature_C__(idx_t)) + 273.15; % average temperature
prop_pma.p = mean(t1.ClassPress_Pa__(idx_t) ./ t1.RefPress_Pa__(idx_t)); % average pressure
sp = get_setpoint(prop_pma,...
    'omega', omega, 'V', V); % get CPMA setpoint information
%===============================================================%



%== Sort based on mass-setpoints ===============================%
m_star = [sp.m_star];
[~,idx_sort] = sort(m_star);

sp = sp(idx_sort);
data = data(:,idx_sort)';
%===============================================================%

end




