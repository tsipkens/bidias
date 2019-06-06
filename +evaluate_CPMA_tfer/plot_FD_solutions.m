

% kk = 301;
kk = 7;

load('cm_inferno.mat');
cm = cm(40:end,:);
figure(1);
contourf(F.r_vec.*100,F.z_vec,F.n_mat{kk},[0:0.05:1.2],'edgecolor','none');
colormap(cm);
caxis([0,1.2]);

axis off;
% set(gca,'Color','w');
% title(['m = ',num2str(m(kk),'%5.2e')]);
% xlabel('r');
% ylabel('z');
% colorbar;






