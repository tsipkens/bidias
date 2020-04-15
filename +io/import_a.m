
% IMPORT_A (EXPERIMENTAL) Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, Arash Naseri, 2019-12-17
%=========================================================================%

function [data,d_star,sp,prop_dma,prop_pma] = import_a(fn_smps,fn_cpma)


%== Read SMPS/CPC file =========================================%
n = io.linecount(fn_smps);
opts = detectImportOptions(fn_smps);
opts.Whitespace = '\b ';
opts.Delimiter = {'\t'};

optsc = opts;
for ii=1:40
    optsc.DataLines = [ii,ii];
    for jj=1:length(optsc.VariableTypes)
        optsc.VariableTypes{jj} = 'string';
    end
    
    ta = readtable(fn_smps,optsc);
    if or(table2array(ta(1,1))=="DMA Inner Radius(cm)",table2array(ta(1,1))=="DMA Inner Radius (cm)");
        nR1 = ii;
    end
    if table2array(ta(1,1))=="Reference Gas Temperature (K)"; nT = ii; end
    if table2array(ta(1,1))=="Date"; nDate = ii; end
end

opts.DataLines = [40,40];
ta = readtable(fn_smps,opts);
ta = table2array(ta);
ta(isnan(ta)) = [];
ncol = length(ta)-1; % number of data samples/scans

opts.DataLines = [nR1,nR1+2];
ta = readtable(fn_smps,opts);
ta = table2array(ta);
prop_dma.R1 = ta(1,2);
prop_dma.R2 = ta(2,2);
prop_dma.L = ta(3,2); % some DMA properties

opts.DataLines = [12,12]; % this is unreliable
ta = readtable(fn_smps,opts);
ta = table2array(ta);
prop_dma.Q_a = ta(1,2)/60/1000;
prop_dma.Q_s = prop_dma.Q_a; % equal flow assumption
prop_dma.Q_c = ta(1,4)/60/1000;
prop_dma.Q_m = prop_dma.Q_c; % equal flow assumption

opts.DataLines = [nT,nT+1];
ta = readtable(fn_smps,opts);
ta = table2array(ta);
prop_dma.T = ta(1,2); % more DMA properties
prop_dma.p = ta(2,2);


optsb = opts;
for ii=1:length(optsb.VariableTypes); optsb.VariableTypes{ii} = 'string'; end
optsb.DataLines = [nDate,nDate+1];
ta = readtable(fn_smps,optsb);
ta = table2array(ta);
for ii=1:ncol
    date_smps(ii,1) = datetime(ta(1,ii+1),'InputFormat','MM/dd/yy')+...
        duration(ta(2,ii+1),'InputFormat','hh:mm:ss');
end


opts.DataLines = [nDate+3,n-26];
t0 = readtable(fn_smps,opts);
t0 = table2array(t0);
data = t0(:,2:(ncol+1));
d_star = t0(:,1);
%==============================================================%



%== Read CPMA file ============================================%
% fid_cpma = fopen(fn_cpma);
% t0 = fgetl(fid_cpma);
% t1 = strsplit(t0,'\t');
% fclose(fid_cpma);

t0 = readtable(fn_cpma);
date_cpma = datetime(t0.Date_Time);

[gr1,gr2] = ndgrid(date_smps,date_cpma);
[~,ind] = max(gr1<gr2,[],2);

omega = t0.ClassSpFB_rad_s__(ind);
V = t0.AdjVoltageFB_V__(ind);

prop_pma = kernel.prop_pma;

for ii=1:length(V)
    prop_pma.T = t0.Temperature_C__(ind(ii))+273.15;
    prop_pma.p = t0.ClassPress_Pa__(ind(ii))/t0.RefPress_Pa__(ind(ii));
    sp(ii) = tfer_pma.get_setpoint(prop_pma,...
        'omega',omega(ii),'V',V(ii));
end

end




