
isom_vec = 10.^(-4:1:3);

for ii=1:length(isom_vec)
    hl = grid_x.overlay_line(...
            [1,log10(6*isom_vec(ii)/(pi*10^3).*1e9)],-3); % line for rho_eff = 0.1
    hl.Color=[1,1,1,0.2];
end
