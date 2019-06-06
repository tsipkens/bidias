

clear;
close all;
% n_case = 21;
% figure(2);
% load('cm_inferno.mat');
% n3 = floor(240/n_case);
% cm_b = flipud(cm(22:n3:256,:));
% set(gca,'ColorOrder',cm_b,'NextPlot','replacechildren');
% m_star_vec = logspace(log10(0.001e-18),log10(10e-18),n_case);
% m_star_vec = 10e-18;
m_star_vec = 0.01e-18;

Rm = 10;

n_rpt = 1;%50;

for ii_m_star = 1:1:length(m_star_vec)

% use m_star = 0.01 to see diffusion, 0.1 to limit diffusion
m_star = m_star_vec(ii_m_star); % mass in kg (1 fg = 1e-18 kg)
% m = logspace(log10(0.1),log10(10),601).*m_star;
m = linspace(0.8,1.2,601).*m_star;
% d = 28e-9.*ones(size(m));
% rho_eff = 6.*m./(pi*d.^3);
rho_eff = 900;%6.*m./(pi*d.^3);
d = (6.*m./(rho_eff.*pi)).^(1/3); % specify mobility diameter


z = 1;

prop = kernel.prop_CPMA('Olfert');
% prop.omega_hat = 1;
D_mod_vec = 1;
prop.D = @(B) prop.D(B)*D_mod_vec(1);

nr = 2000;
dr = (prop.r2-prop.r1)/nr;
r_vec = (prop.r1+dr/2):dr:(prop.r2-dr/2);

for tt=1:n_rpt
tic;
[tfer_FD,~,F] = kernel.tfer_CPMA_FD(m_star,...
    m,d,1,prop,'Rm',Rm);
t(1,tt) = toc;
end

%-------------------------------------------------------------------------%
%-- Transfer functions for different cases -------------------------------%
%-- Setup for centriputal force ------------------------------------------%
if ~exist('d','var')
    B = kernel.mp2zp(m,z,prop.T,prop.p);
else
    B = kernel.dm2zp(d,z,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B).*z;
sig = sqrt(2.*prop.L.*D./prop.v_bar);
D0 = D.*prop.L/(prop.del^2*prop.v_bar); % dimensionless diffusion coeff.



%-- Particle tracking approaches -----------------------------------------%
for tt=1:n_rpt

%-- Plug flow ------------------------------------------------------------%
%-- Method A ------------------------------%
tic;
[tfer_A,G0_A] = kernel.tfer_CPMA_A(m_star,m,d,1,prop,'Rm',Rm);
t(2,tt) = toc;

%-- Method A, Ehara et al. ----------------%
tfer_A_Ehara = kernel.tfer_CPMA_A_Ehara(m_star,m,d,1,prop,'Rm',Rm);

%-- Method B ------------------------------%
tic;
[tfer_B,G0_B] = kernel.tfer_CPMA_B(m_star,m,d,1,prop,'Rm',Rm);
t(3,tt) = toc;

[tfer_B2] = kernel.tfer_CPMA_B(m_star,m,d,2,prop,'Rm',Rm);

%-- Method C ------------------------------%
tic;
[tfer_C,G0_C] = kernel.tfer_CPMA_C(m_star,m,d,1,prop,'Rm',Rm);
t(4,tt) = toc;

%-- Method D ------------------------------%
tic;
[tfer_D,G0_D] = kernel.tfer_CPMA_D(m_star,m,d,1,prop,'Rm',Rm);
t(5,tt) = toc;

%-- Method E ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_E,G0_E] = kernel.tfer_CPMA_E(m_star,m,d,z,prop,'Rm',Rm);
    t(6,tt) = toc;
end

%-- Method F ------------------------------%
tic;
[tfer_F,G0_F] = kernel.tfer_CPMA_F(m_star,m,d,z,prop,'Rm',Rm);
t(7,tt) = toc;


%-- Parabolic flow -------------------------------------------------------%
%-- Method A ------------------------------%
tic;
[tfer_A_pb,G0_A_pb] = kernel.tfer_CPMA_A_pb(m_star,m,d,z,prop,'Rm',Rm);
t(8,tt) = toc;

%-- Method B ------------------------------%
tic;
[tfer_B_pb,G0_B_pb] = kernel.tfer_CPMA_B_pb(m_star,m,d,z,prop,'Rm',Rm);
t(9,tt) = toc;

%-- Method E ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_E_pb,G0_E_pb] = kernel.tfer_CPMA_E_pb(m_star,m,d,z,prop,'Rm',Rm);
    t(10,tt) = toc;
end


%-- Diffusive transfer functions -----------------------------------------%
%-- Method A ------------------------------%
tic;
tfer_A_diff = kernel.tfer_CPMA_A_diff(m_star,m,d,z,prop,'Rm',Rm);
t(11,tt) = toc;

%-- Method B -------------------------------%
tic;
tfer_B_diff = kernel.tfer_CPMA_B_diff(m_star,m,d,z,prop,'Rm',Rm);
t(12,tt) = toc;

%-- Method C -----------------------------%
tic;
tfer_C_diff = kernel.tfer_CPMA_C_diff(m_star,m,d,z,prop,'Rm',Rm);
t(13,tt) = toc;

%-- Method D --------------------------------%
tic;
tfer_D_diff = kernel.tfer_CPMA_D_diff(m_star,m,d,z,prop,'Rm',Rm);
t(14,tt) = toc;

%-- Method E --------------------------------%
if prop.omega_hat==1
    tic;
    tfer_E_diff = kernel.tfer_CPMA_E_diff(m_star,m,d,z,prop,'Rm',Rm);
    t(15,tt) = toc;
end


%-- Triangle approx. -----------------------%
tic;
tfer_tri = kernel.tfer_CPMA_tri(m_star,m,d,z,prop,'Rm',Rm);
t(18,tt) = toc;

end



%%

figure(2);
plot(m,tfer_A);
% hold on;
% plot(m,tfer_A_Ehara);
% plot(m,tfer_A_pb);
% plot(m,tfer_A_diff);
hold on;
plot(m,tfer_B);
% hold on;
plot(m,tfer_B_diff);
% plot(m./m_star(1),tfer_B_diff);
% hold on;
% plot(m,tfer_B_pb);
% plot(m,tfer_B_pb_alt);
% plot(m,tfer_C);
% plot(m,tfer_C_diff);
% hold on;
% plot(m,tfer_D_diff);
% plot(m,tfer_D);
% hold on;
% plot(m,tfer_E,'r');
% hold on;
% plot(m,tfer_E_diff,'r');
% plot(m,tfer_E_pb);
% plot(m,tfer_F);
plot(m,tfer_tri);
% hold on;
plot(m,min(tfer_FD,1),'k');
% plot(m./m_star(1),min(tfer_FD,1));
hold off;

% xlim([0.85,1.15]);
% xlim([0.85,1.15].*m_star(1));
ylim([0,1.2]);
% ylim([0,0.8]);


end

%{
figure(3);
plot(m,tfer_A-min(tfer_FD,1));
hold on;
plot(m,tfer_B-min(tfer_FD,1));
plot(m,tfer_C-min(tfer_FD,1));
plot(m,tfer_D-min(tfer_FD,1));
if prop.omega_hat==1; plot(m,tfer_E-min(tfer_FD,1)); end
plot(m,tfer_F-min(tfer_FD,1));
plot(m,0.*tfer_FD);
hold off;
xlim([0.85,1.15].*m_star(1));
%}

%{
figure(3);
plot(m,tfer_A_pb-min(tfer_FD,1));
hold on;
plot(m,tfer_B_pb-min(tfer_FD,1));
if prop.omega_hat==1; plot(m,tfer_E_pb-min(tfer_FD,1)); end
plot(m,0.*tfer_FD);
hold off;
xlim([0.85,1.15].*m_star(1));
%}

%-{
figure(3);
plot(m,tfer_A_diff-min(tfer_FD,1));
hold on;
plot(m,tfer_B_diff-min(tfer_FD,1));
plot(m,tfer_C_diff-min(tfer_FD,1));
plot(m,tfer_D_diff-min(tfer_FD,1));
% plot(m,tfer_E_diff-min(tfer_FD,1));
% plot(m,tfer_F-min(tfer_FD,1));
plot(m,0.*tfer_FD);
hold off;
xlim([0.85,1.15].*m_star(1));
%}


% bar(log10(mean(t,2))+5);


%%
%{
vec = 290:454;
chi.A = norm(tfer_A(vec)-tfer_FD(vec));
chi.B = norm(tfer_B(vec)-tfer_FD(vec));
chi.C = norm(tfer_C(vec)-tfer_FD(vec));
chi.D = norm(tfer_D(vec)-tfer_FD(vec));
% chi.E = norm(tfer_E(vec)-tfer_FD(vec));
chi.F = norm(tfer_F(vec)-tfer_FD(vec));
chi.A_pb = norm(tfer_A_pb(vec)-tfer_FD(vec));
chi.B_pb = norm(tfer_B_pb(vec)-tfer_FD(vec));
% chi.E_pb = norm(tfer_E_pb(vec)-tfer_FD(vec));
chi.A_diff = norm(tfer_A_diff(vec)-tfer_FD(vec));
chi.B_diff = norm(tfer_B_diff(vec)-tfer_FD(vec));
chi.C_diff = norm(tfer_C_diff(vec)-tfer_FD(vec));
chi.D_diff = norm(tfer_D_diff(vec)-tfer_FD(vec));
% chi.E_diff = norm(tfer_E_diff(vec)-tfer_FD(vec));
chi.tri = norm(tfer_tri-tfer_FD);

figure(4);
chi_names = fieldnames(chi);
chi_vals = zeros(length(chi_names),1);
for ii=1:length(chi_names)
    chi_vals(ii) = chi.(chi_names{ii});
end

bar(chi_vals);
set(gca,'xticklabel',chi_names);
%}


figure(2);
