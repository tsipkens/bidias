
% IMPORT_A Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, Arash Naseri, 2019-12-17
%=========================================================================%

function [sp,prop_dma,prop_pma,out] = import_a(fn_smps,fn_cpma)

fid_smps = fopen(fn_smps);
fid_cpma = fopen(fn_cpma);

fclose(fid_smps);
fclose(fid_cpma);

sp.d_star = [];
sp.m_star = [];

end


