

clear;

n_rpt = 1;%50
t = [];

d_star = [];
m_star = 1e-18; % mass in kg (1 fg = 1e-18 kg)
% m = [0.85,0.9,0.95,1,1.05,1.1,1.15].*m_star;
m = logspace(log10(0.5),log10(1.5),599).*m_star;
z = 1;

prop = kernel.prop_CPMA('Olfert');
prop.omega_hat = 1;
D_mod_vec = logspace(log10(0.01),log10(10),19);
% prop.D = prop.D*D_mod_vec(19);
prop.Omega = prop.D*prop.L/(prop.del^2*prop.v_bar); % dimensionless diffusion coeff.

nr = 2000;
dr = (prop.r2-prop.r1)/nr;
r_vec = (prop.r1+dr/2):dr:(prop.r2-dr/2);


for tt=1:n_rpt
tic;
tfer_FD = kernel.tfer_CPMA_FD(m_star,...
    m,[],1,prop)';
% tfer_FD_RK = kernel.tfer_CPMA_FD_RK(m_star,...
%     m,[],1,prop)';
t(1,tt) = toc;
end

%-------------------------------------------------------------------------%
%-- Transfer functions for different cases -------------------------------%
kernel.prop_CPMA_plus;
sig = 2*sqrt(prop.L*prop.D/prop.v_bar);


%-- Particle tracking approaches -----------------------------------------%
for tt=1:n_rpt

%-- Plug flow ------------------------------------------------------------%
%-- Method A ------------------------------%
tic;
[tfer_A,G0_A] = kernel.tfer_CPMA_A(prop,lam,rs);
t(2,tt) = toc;

%-- Method A, Ehara et al. ----------------%
rho_s = (rs-prop.rc)/prop.del;
rho_1 = (prop.r1-prop.rc)/prop.del;
rho_2 = (prop.r2-prop.rc)/prop.del;
tfer_A_Ehara = ((1-rho_s)+(1+rho_s).*exp(-lam))./2.*and(1<rho_s,rho_s<coth(lam./2))+...
    exp(-lam).*and(-1<rho_s,rho_s<1)+...
    ((1+rho_s)+(1-rho_s).*exp(-lam))./2.*and(-coth(lam./2)<rho_s,rho_s<-1);

%-- Method A, alternate --------------------%
% tic;
% f_min = @(rL,r0,ii) exp(lam(ii)).*(r0-rs(ii))-(rL-rs(ii));
% condit = lam<0;
% G0_A = @(r) kernel.G_fun(f_min,r,rs,prop.r1,prop.r2,condit);
% G1_A = G0_A(prop.r1);
% G2_A = G0_A(prop.r2);
% ra = min(prop.r2,max(prop.r1,G1_A));
% rb = min(prop.r2,max(prop.r1,G2_A));
% tfer_A_alt = (1/(2*prop.del)).*(rb-ra);
% t(9,tt) = toc;

%-- Method B ------------------------------%
tic;
[tfer_B,G0_B] = kernel.tfer_CPMA_B(prop,C3,C4);
t(3,tt) = toc;

%-- Method C ------------------------------%
tic;
[tfer_C,G0_C] = kernel.tfer_CPMA_C(prop,lam,rs,C1,C2);
t(4,tt) = toc;

%-- Method D ------------------------------%
tic;
[tfer_D,G0_D] = kernel.tfer_CPMA_D(prop,C3,C4,C5);
t(5,tt) = toc;

%-- Method E ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_E,G0_E] = kernel.tfer_CPMA_E(prop,lam,C0,m);
    t(6,tt) = toc;
end

%-- Method F ------------------------------%
tic;
[tfer_F,G0_F] = kernel.tfer_CPMA_F(prop,tau,C0,m,rs);
t(7,tt) = toc;



%-- Parabolic flow -------------------------------------------------------%
%-- Method A ------------------------------%
tic;
[tfer_A_pb,G0_A_pb] = kernel.tfer_CPMA_A_pb(prop,lam,rs);
t(8,tt) = toc;

%-- Method B ------------------------------%
tic;
[tfer_B_pb,G0_B_pb] = kernel.tfer_CPMA_B_pb(prop,C3,C4,rs);
t(9,tt) = toc;

%-- Method E ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_E_pb,G0_E_pb] = kernel.tfer_CPMA_E_pb(prop,tau,C0,m,rs);
    t(10,tt) = toc;
end



%-- Diffusive transfer functions -----------------------------------------%
rho_fun = @(G,r) (G-r)./(sqrt(2)*sig);
K_fun = @(G,r) ...
    (G-r).*erf(rho_fun(G,r))+...
    sig*sqrt(2/pi).*exp(-rho_fun(G,r).^2);

K1_fun = @(G,r) ((G-r)+...
    1/(prop.del^2).*(-prop.rc^2.*(G-r)+prop.rc.*(G.^2-r.^2-sig^2)+r.*sig^2-(G.^3-r.^3)./3)).*...
    erf(rho_fun(G,r));
K2_fun = @(G,r) sig*sqrt(2/pi).*...
    (1+1/(prop.del).*(-prop.rc^2+prop.rc.*(G+r)-(G+r).^2./3-2*sig^2/3)).*...
    exp(-rho_fun(G,r).^2);
K_pb_fun = @(G,r) K1_fun(G,r)+K2_fun(G,r);

%-- Method A ------------------------------%
tic;
K22 = K_fun(G0_A(prop.r2),prop.r2);
K21 = K_fun(G0_A(prop.r2),prop.r1);
K12 = K_fun(G0_A(prop.r1),prop.r2);
K11 = K_fun(G0_A(prop.r1),prop.r1);
tfer_A_diff = -1/(4*prop.del).*(K22-K12-K21+K11);
t(11,tt) = toc;

%-- Method B -------------------------------%
tic;
K22 = K_fun(G0_B(prop.r2),prop.r2);
K21 = K_fun(G0_B(prop.r2),prop.r1);
K12 = K_fun(G0_B(prop.r1),prop.r2);
K11 = K_fun(G0_B(prop.r1),prop.r1);
tfer_B_diff = -1/(4*prop.del).*(K22-K12-K21+K11);
t(12,tt) = toc;

%-- Method C -----------------------------%
tic;
K22 = K_fun(G0_C(prop.r2),prop.r2);
K21 = K_fun(G0_C(prop.r2),prop.r1);
K12 = K_fun(G0_C(prop.r1),prop.r2);
K11 = K_fun(G0_C(prop.r1),prop.r1);
tfer_C_diff = -1/(4*prop.del).*(K22-K12-K21+K11);
t(13,tt) = toc;

%-- Method D --------------------------------%
tic;
K22 = K_fun(real(G0_D(prop.r2)),prop.r2);
K21 = K_fun(real(G0_D(prop.r2)),prop.r1);
K12 = K_fun(real(G0_D(prop.r1)),prop.r2);
K11 = K_fun(real(G0_D(prop.r1)),prop.r1);
tfer_D_diff = -1/(4*prop.del).*(K22-K12-K21+K11);
t(14,tt) = toc;

%-- Method E --------------------------------%
if prop.omega_hat==1
    tic;
    K22 = K_fun(G0_E(prop.r2),prop.r2);
    K21 = K_fun(G0_E(prop.r2),prop.r1);
    K12 = K_fun(G0_E(prop.r1),prop.r2);
    K11 = K_fun(G0_E(prop.r1),prop.r1);
    tfer_E_diff = -1/(4*prop.del).*(K22-K12-K21+K11);
    t(15,tt) = toc;
end

%-- Method A, parabolic ---------------------%
% tic;
% K22 = K_fun(G0_A_pb(prop.r2),prop.r2);
% K21 = K_fun(G0_A_pb(prop.r2),prop.r1);
% K12 = K_fun(G0_A_pb(prop.r1),prop.r2);
% K11 = K_fun(G0_A_pb(prop.r1),prop.r1);
% tfer_A_pb_diff = -3/(8*prop.del).*(K22-K12-K21+K11);
% t(16,tt) = toc;

%-- Method B, parabolic --------------------%
% tic;
% K22 = K_fun(G0_B_pb(prop.r2),prop.r2);
% K21 = K_fun(G0_B_pb(prop.r2),prop.r1);
% K12 = K_fun(G0_B_pb(prop.r1),prop.r2);
% K11 = K_fun(G0_B_pb(prop.r1),prop.r1);
% tfer_B_pb_diff = -3/(8*prop.del).*(K22-K12-K21+K11);
% t(17,tt) = toc;



%-- Triangle approx. -----------------------%
L_max = 1;
md = 0.1;
mc = 1;
tfer_tri_fun = @(x) x(1).*((m>(x(3).*1e-18-x(2).*1e-18)).*(m<=x(3).*1e-18).*(m-(x(3).*1e-18-x(2).*1e-18))./(x(2).*1e-18)+...
    (m>x(3).*1e-18).*(m<(x(3).*1e-18+x(2).*1e-18)).*(1+(x(3).*1e-18-m)./(x(2).*1e-18)));
x1 = lsqnonlin(@(x) tfer_tri_fun(x)-tfer_FD,[L_max,md,mc],...
    [0,0,0],[1,inf,inf]);

tic;
tfer_tri = tfer_tri_fun(x1);
t(18,tt) = toc;

end

% figure(2);
% plot(tfer_A);
% hold on;
% plot(tfer_A_diff);
% plot(tfer_A_alt,'--');
% plot(tfer_A_diff_alt,'--');
% hold off;



%%
load('cm_inferno.mat');
cm = cm(22:13:256,:);
set(gca,'colororder',cm);

figure(2);
% plot(m,tfer_A);
% hold on;
% plot(m,tfer_A_Ehara);
% plot(m,tfer_A_pb);
% plot(m,tfer_A_diff);
% hold on;
% plot(m,tfer_A_pb_diff)
plot(m,tfer_B);
hold on;
plot(m,tfer_B_diff);
% hold on;
plot(m,tfer_B_pb);
% plot(m,tfer_B_pb_diff)
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
plot(m,min(tfer_FD,1));
hold off;

xlim([0.85,1.15].*1e-18);
ylim([0,1.2]);
% ylim([0,0.8]);

% plot(r_vec,(2.*C6.*prop.v_bar).*t1);
% hold on;
% plot(r_vec,(2.*C6.*prop.v_bar).*t0,'--')
% hold off;

%{
figure(3);
plot(m,tfer_A_diff-min(tfer_FD,1));
hold on;
plot(m,tfer_B_diff-min(tfer_FD,1));
plot(m,tfer_C_diff-min(tfer_FD,1));
plot(m,tfer_D_diff-min(tfer_FD,1));
plot(m,tfer_E_diff-min(tfer_FD,1));
% plot(m,tfer_F-min(tfer_FD,1));
plot(m,0.*tfer_FD);
hold off;
xlim([0.85,1.15].*1e-18);
%}


% bar(log10(mean(t,2))+5);


%%
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



