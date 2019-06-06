


figure(1);

d_star = 20*1e-9;
d_plot = ((0.8*d_star):0.01*1e-9:(1.3*d_star));
d = d_plot;

% olfert.DMA_D_Stolzenburg

Z1 = kernel.dm2zp(20e-9,1,298,1);
Z2 = olfert.dp2elecmob(20e-9,1,1,298);

% f_DMA_old = kernel.f_DMA_old(d_star,d_plot,1);
% plot(kernel.dp2zp(d_plot)./kernel.dp2zp(d_star),f_DMA_old);

% [f,Zp_tilde] = kernel.f_DMA(d_star,d_plot,1,1,'Buckley');
% hold on;
% plot(Zp_tilde,f);
% hold off;
% 
% [f,Zp_tilde] = kernel.f_DMA(d_star,d_plot,1,1,'plug');
% hold on;
% plot(Zp_tilde,f);
% hold off;

hold on;
[f,Zp_tilde] = kernel.tfer_DMA(d_star,d_plot,1,1,'Olfert');
hold on;
plot(Zp_tilde,f);
hold off;

[f,Zp_tilde] = kernel.tfer_DMA(d_star,d_plot,1,0);
hold on;
plot(Zp_tilde,f,'-');
hold off;

xlim([0.6,1.4]);

