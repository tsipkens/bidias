
% IMPORT_A (EXPERIMENTAL) Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, Arash Naseri, 2019-12-17
%=========================================================================%

function [data,d_star,sp,prop_dma,prop_pma] = import_a(fn_smps,fn_cpma)


%== Read SMPS/CPC file =========================================%
[data,d_star,prop_dma,time_dma] = import_dma(fn_smps);
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




