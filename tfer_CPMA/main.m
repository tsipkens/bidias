

clear;
close all;

Rm = 10; % resolution of transfer functions

m_star = 0.01e-18; % mass in kg (1 fg = 1e-18 kg)
m = linspace(0.8,1.2,601).*m_star; % vector of mass

z = 1; % integer charge state

rho_eff = 900; % effective density
d = (6.*m./(rho_eff.*pi)).^(1/3);
    % specify mobility diameter vector with constant effective density

prop = kernel_CPMA.prop_CPMA('Olfert'); % get properties of the CPMA
% prop.omega_hat = 1; % modify the rotational speed ratio



%%
%-------------------------------------------------------------------------%
%-- Finite difference solutions ------------------------------------------%
tic;
[tfer_FD,~,F] = kernel_CPMA.tfer_CPMA_FD(m_star,...
    m,d,1,prop,'Rm',Rm);
t(1) = toc;



%%
%-------------------------------------------------------------------------%
%-- Transfer functions for different cases -------------------------------%
%-- Setup for centriputal force ------------------------------------------%
if ~exist('d','var')
    B = kernel_CPMA.mp2zp(m,z,prop.T,prop.p);
else
    B = kernel_CPMA.dm2zp(d,z,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B).*z;
sig = sqrt(2.*prop.L.*D./prop.v_bar);
D0 = D.*prop.L/(prop.del^2*prop.v_bar); % dimensionless diffusion coeff.


%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method A ------------------------------%
tic;
[tfer_A,G0_A] = kernel_CPMA.tfer_CPMA_A(m_star,m,d,1,prop,'Rm',Rm);
t(2) = toc;

%-- Method A, Ehara et al. ----------------%
tfer_A_Ehara = kernel_CPMA.tfer_CPMA_A_Ehara(m_star,m,d,1,prop,'Rm',Rm);

%-- Method B ------------------------------%
tic;
[tfer_B,G0_B] = kernel_CPMA.tfer_CPMA_B(m_star,m,d,1,prop,'Rm',Rm);
t(3) = toc;

[tfer_B2] = kernel_CPMA.tfer_CPMA_B(m_star,m,d,2,prop,'Rm',Rm);

%-- Method C ------------------------------%
tic;
[tfer_C,G0_C] = kernel_CPMA.tfer_CPMA_C(m_star,m,d,1,prop,'Rm',Rm);
t(4) = toc;

%-- Method D ------------------------------%
tic;
[tfer_D,G0_D] = kernel_CPMA.tfer_CPMA_D(m_star,m,d,1,prop,'Rm',Rm);
t(5) = toc;

%-- Method E ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_E,G0_E] = kernel_CPMA.tfer_CPMA_E(m_star,m,d,z,prop,'Rm',Rm);
    t(6) = toc;
end

%-- Method F ------------------------------%
tic;
[tfer_F,G0_F] = kernel_CPMA.tfer_CPMA_F(m_star,m,d,z,prop,'Rm',Rm);
t(7) = toc;


%-- Parabolic flow -------------------------------------------------------%
%-- Method A ------------------------------%
tic;
[tfer_A_pb,G0_A_pb] = kernel_CPMA.tfer_CPMA_A_pb(m_star,m,d,z,prop,'Rm',Rm);
t(8) = toc;

%-- Method B ------------------------------%
tic;
[tfer_B_pb,G0_B_pb] = kernel_CPMA.tfer_CPMA_B_pb(m_star,m,d,z,prop,'Rm',Rm);
t(9) = toc;

%-- Method E ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_E_pb,G0_E_pb] = kernel_CPMA.tfer_CPMA_E_pb(m_star,m,d,z,prop,'Rm',Rm);
    t(10) = toc;
end


%-- Diffusive transfer functions -----------------------------------------%
%-- Method A ------------------------------%
tic;
tfer_A_diff = kernel_CPMA.tfer_CPMA_A_diff(m_star,m,d,z,prop,'Rm',Rm);
t(11) = toc;

%-- Method B -------------------------------%
tic;
tfer_B_diff = kernel_CPMA.tfer_CPMA_B_diff(m_star,m,d,z,prop,'Rm',Rm);
t(12) = toc;

%-- Method C -------------------------------%
tic;
tfer_C_diff = kernel_CPMA.tfer_CPMA_C_diff(m_star,m,d,z,prop,'Rm',Rm);
t(13) = toc;

%-- Method D --------------------------------%
tic;
tfer_D_diff = kernel_CPMA.tfer_CPMA_D_diff(m_star,m,d,z,prop,'Rm',Rm);
t(14) = toc;

%-- Method E --------------------------------%
if prop.omega_hat==1
    tic;
    tfer_E_diff = kernel_CPMA.tfer_CPMA_E_diff(m_star,m,d,z,prop,'Rm',Rm);
    t(15) = toc;
end

%-- Method F --------------------------------%
tic;
tfer_F_diff = kernel_CPMA.tfer_CPMA_F_diff(m_star,m,d,z,prop,'Rm',Rm);
t(16) = toc;

%-- Triangle approx. -----------------------%
tic;
tfer_tri = kernel_CPMA.tfer_CPMA_tri(m_star,m,d,z,prop,'Rm',Rm);
t(17) = toc;



%%
m_plot = m./m_star;

figure(2);
plot(m_plot,tfer_A);
hold on;
% plot(m_plot,tfer_A_Ehara);
% plot(m_plot,tfer_A_diff);
% plot(m_plot,tfer_A_pb);
plot(m_plot,tfer_B);
plot(m_plot,tfer_B_diff);
% plot(m_plot,tfer_B_pb);
% plot(m_plot,tfer_C);
% plot(m_plot,tfer_C_diff);
% plot(m_plot,tfer_D);
% plot(m_plot,tfer_D_diff);
% plot(m_plot,tfer_E,'r');
% plot(m_plot,tfer_E_diff,'r');
% plot(m_plot,tfer_E_pb);
% plot(m_plot,tfer_F);
% plot(m_plot,tfer_F_diff);
% plot(m_plot,tfer_tri);
plot(m_plot,min(tfer_FD,1),'k');
hold off;

ylim([0,1.2]);
xlim(1+[-1.5/Rm,1.5/Rm]);




