

olfert.test_file;
tfer_Olfert_Df = trans_Df;
tfer_Olfert_nDf = trans;
tfer_Olfert_tri = Omega;

kernel.sim_CPMA_v4;
tfer_int1 = tfer;
%%
prop_CPMA = kernel.prop_CPMA('Olfert');
tfer_int1 = kernel.tfer_CPMA(mc,m,[],1,prop_CPMA);
prop_CPMA.omega_hat = 1;
tfer_int2 = kernel.tfer_CPMA(mc,m,[],1,prop_CPMA);
tfer_int_APM2 = kernel.tfer_APM(mc,m,[],1,prop_CPMA);

prop_APM = kernel.prop_CPMA('Buckley');
tfer_int_APM3 = kernel.tfer_APM(mc,m,[],1,prop_APM);
tfer_int3 = kernel.tfer_CPMA(mc,m,[],1,prop_APM);

figure(1);
plot(m,tfer_Olfert_Df);
hold on;
plot(m,tfer_Olfert_nDf);
plot(m,tfer_Olfert_tri);
plot(m,tfer_int1,'--');

plot(m,tfer_int2,'r.');
plot(m,tfer_int_APM2,'r');
plot(m,tfer_int_APM3,'k');
plot(m,tfer_int3,'.k');
hold off;

%%

prop = kernel.prop_CPMA('Olfert');
m_star_vec = logspace(log10(0.1),log10(100),19).*1e-18;
m = logspace(log10(0.01),log10(1000),600).*1e-18;

figure(2);
load('..\Program - LII\LII Program 3.9\+CENIDE\viridis.mat');
cm = cm(1:floor(256/18):end,:);
tfer0 = [];
h = [];
clf;
semilogx(mean(m_star_vec),0.1);
for ii=1:length(m_star_vec)
    tfer0{ii} = kernel.tfer_CPMA(m_star_vec(ii),m,[],1,prop);
    hold on;
    h(ii) = semilogx(m,tfer0{ii});
    hold off;
    
    if m_star_vec(ii)==1e-17
        set(h(ii),'color',[0,0,0]);
    else
        set(h(ii),'color',cm(ii,:));
    end
end
% hold on;
% plot([0.8e-17,0.8e-17],[0,0.9]);
% plot([0.9e-17,0.9e-17],[0,0.9]);
% plot([1.0e-17,1.0e-17],[0,0.9]);
% plot([1.1e-17,1.1e-17],[0,0.9]);
% plot([1.2e-17,1.2e-17],[0,0.9]);
% hold off;
xlim([5e-20,5e-16]);

%%
t0 = kernel.tfer_CPMA(1e-17,1.2.*1e-17,[],1,prop);
caxis([0,1.1]);
colorbar;



