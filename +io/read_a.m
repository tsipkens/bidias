
% READ_A  Read tandem CPMA-SMPS data using SMPS data export and a CPMA log file.
% Author: Timothy Sipkens, Arash Naseri, 2019-12-17
%=========================================================================%

function [sp,prop_dma,prop_pma,out] = read_a(fn_smps,fn_cpma)

fid_smps = fopen(fn_smps);
fid_cpma = fopen(fn_cpma);

fclose(fid_smps);
fclose(fid_cpma);

end


