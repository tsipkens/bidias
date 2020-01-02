
% IMPORT_A Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, Arash Naseri, 2019-12-17
%=========================================================================%

function [data,d_star,sp,prop_dma,prop_pma] = import_a(fn_smps,fn_cpma)

fid_smps = fopen(fn_smps);

t0 = textscan(fid_smps,'%s','delimiter','\t');
prop_dma = struct(); % properties structure require for analysis
prop_dma_alt = struct(); % alternate properties read directly from file

n = length(t0{1});


%== Read SMPS file ============================================%
for ii=1:2:n % read DMA properties
    t1 = t0{1}{ii};
    if strcmp(t0{1}{ii},'Sample #'); break; end
        % finished reading DMA properties
    
    t1 = strrep(t1,' ',''); % remove invalid field name characters
    t1 = strrep(t1,'(','_');
    t1 = strrep(t1,')','');
    t1 = strrep(t1,'*','_');
    t1 = strrep(t1,'/','_per_');
    t1 = strrep(t1,'#','No');
    
    prop_dma_alt.(t1) = t0{1}{ii+1};
end


%-- Determine the number of samples in the file ----%
for jj=ii+1:1:n
    if strcmp(t0{1}{jj},'Date'); break; end
end
n_samples = jj-ii-1;


%-- Read time and data -----------------------------%
for kk=1:n_samples
    smps_time(kk) = datetime(...
        [t0{1}{kk+jj},' ',...
        t0{1}{kk+jj+n_samples+1}],...
        'InputFormat','MM/dd/yy HH:mm:ss',...
        'PivotYear',2000);
end

%-- Read CPC data ---------------------------------------------%
kk = 0;
for jj=(jj+2*(n_samples+1)+1):(n_samples+1):n
    if strcmp(t0{1}{jj},'Scan Up Time(s)'); break; end
        % finished reading data
    
    kk = kk+1;
    d_star(kk) = str2double(t0{1}{jj});
    
    for ll=1:n_samples
        data(kk,ll) = str2double(t0{1}{jj+ll});
    end
end

%-- Read in flow (assume the same over all runs in file) ------%
prop_dma.Q_c = str2double(t0{1}{jj+5*(n_samples+1)+1})/60/1000; % sheath flow [m^3/s]
prop_dma.Q_a = str2double(t0{1}{jj+6*(n_samples+1)+1})/60/1000; % aerosol flow [m^3/s]
prop_dma.Q_s = prop_dma.Q_a; % sample flow [m^3/s] (assume equal flow)
prop_dma.Q_m = prop_dma.Q_c; % exhaust flow [m^3/s] (assume equal flow)

%-- Reassign properties in alternate struct -------------------%
prop_dma.R1 = str2double(prop_dma_alt.DMAInnerRadius_cm);
prop_dma.R2 = str2double(prop_dma_alt.DMAOuterRadius_cm);
prop_dma.T = str2double(prop_dma_alt.ReferenceGasTemperature_K);
prop_dma.p = str2double(prop_dma_alt.ReferenceGasPressure_kPa)/101.325;

fclose(fid_smps);



%== Read CPMA file ============================================%
% fid_cpma = fopen(fn_cpma);
% t0 = fgetl(fid_cpma);
% t1 = strsplit(t0,'\t');
% fclose(fid_cpma);

t0 = readtable(fn_cpma);
t1 = t0.Date_Time;

[gr1,gr2] = ndgrid(smps_time,t1);
[~,ind] = max(gr1<gr2,[],2);

omega = t0.ClassSpFB_rad_s__(ind);
V = t0.AdjVoltageFB_V__(ind);

prop_pma = kernel.prop_pma('fn18');
prop_pma.mass_mob_pref = 0.0612; % assume CPMA uses soot properties
prop_pma.mass_mob_exp = 2.48;

for ii=1:length(V)
    prop_pma.T = t0.Temperature_C__(ind(ii))+273.15;
    prop_pma.p = t0.ClassPress_Pa__(ind(ii))/t0.RefPress_Pa__(ind(ii));
    sp(ii) = tfer_pma.get_setpoint(prop_pma,...
        'omega',omega(ii),'V',V(ii));
end

end


