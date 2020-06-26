
% IMPORT_2  Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, 2019-12-17
%=========================================================================%

function [data,d_star,sp,prop_dma,prop_pma] = import_2(fn_smps,fn_cpma,prop_pma)


%== Read SMPS/CPC file =========================================%
[data,d_star,prop_dma,time_dma] = io.import_dma(fn_smps);
%===============================================================%



%== Read CPMA file =============================================%
addpath('tfer_pma');

opts = detectImportOptions(fn_cpma); % default options
t0 = readtable(fn_cpma, opts); % read CPMA file
time_cpma = datetime(t0.Date_Time); % get times to compare against DMA data

[gr1,gr2] = ndgrid(time_dma, time_cpma); % grid of CPMA/SMPS times
[~,ind] = max(gr1<gr2,[],2); % get first CPMA time corresponding to DMA scan

if any(ind==1)
    warning('CPMA data files do not contain correct times. Only DMA data output');
    sp = struct();
    return;
end

omega = t0.ClassSpFB_rad_s__(ind); % speed of CPMA (no averaging)
V = t0.AdjVoltageFB_V__(ind); % voltage of CPMA (no averaging)

% alternate voltages / speeds (unused)
V_raw = t0.RawVoltageFB_V__(ind);
omega1 = t0.InnerSpFB_rad_s__(ind);
omega2 = t0.OuterSpFB_rad_s__(ind);

% if prop_pma is not specified
if ~exist('prop_pma','var'); prop_pma = []; end
if isempty(prop_pma); prop_pma = kernel.prop_pma; end


prop_pma.T = mean(t0.Temperature_C__(ind)) + 273.15; % average temperature
prop_pma.p = mean(t0.ClassPress_Pa__(ind) ./ t0.RefPress_Pa__(ind)); % average pressure
sp = get_setpoint(prop_pma,...
    'omega', omega, 'V', V); % get CPMA setpoint information
%===============================================================%


end




