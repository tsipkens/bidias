
isorho_vec = 10.^(-1:1:7);

for ii=1:length(isorho_vec)
    hl = grid_x.overlay_line(...
            [1,log10(isorho_vec(ii)*10^3./1e9)],3); % line for rho_eff = 0.1
    hl.Color=[1,1,1,0.2];
end
