
function [Lambda,prop] = tfer_CPMA(m_star,m,d,z,prop,opts)
% TFER_CPMA Evaluates the transfer function for a centrifugal particle mass analyzer (CPMA).
% Author:       Timothy Sipkens, 2018-12-27
% Based on:     Buckley et al., J. Aerosol Sci. (2008) and Olfert group
% 
%--------------------------------------------------------------------------%
% Inputs:
%   m_star      Mass corresponding to the measurement set point of the APM
%   d           Particle mobility diameter, can be vector [nm]
%   m           Particle mass, can be vector of same length as d
%   z           Integer charge state, scalar
%   prop        CPMA device settings (optional)
%   opts        Options structure (optional)
%
% Outputs:
%   Lambda      CPMA transfer function
%   prop        CPMA device settings
%-------------------------------------------------------------------------%

if ~exist('prop','var') % if order not specified
    prop = kernel.prop_CPMA;
elseif isempty(prop)
    prop = kernel.prop_CPMA;
end

if ~exist('opts','var')
    opts = [];
end


%-- Discretize the space -------------------------------------------------%
nr = 80;
dr = (prop.r2-prop.r1)/nr;
r_vec = (prop.r1+dr/2):dr:(prop.r2-dr/2);

nz = 80;
dz = prop.L/(nz-1);
% z_vec = 0:dz:prop.L; % vector of z locations, used for plotting only


%-- Set up physical parameters -------------------------------------------%
e = 1.6e-19; % electron charge [C]
q = z.*e;
kB = 1.3806488e-23; % Boltzmann's constant


%-- Parameters related to CPMA geometry ----------------------------------%
del = (prop.r2-prop.r1)/2;
rc = (prop.r2+prop.r1)/2;
r_hat = prop.r1/prop.r2;


%-- Evaluate set point parameters ----------------------------------------%
if ~isfield(prop,'omega') % if rotational speed is not explicitly specified, use formulas from Olfert lab
    n_B = -0.6436;
    B_star = kernel.mp2zp(m_star,1,prop.T,prop.p); % use z = 1 for CPMA setpoint
    Rm = 10;
    m_max = m_star*(1/Rm+1);
    prop.omega = sqrt(prop.Q/(m_star*B_star*2*pi*rc^2*prop.L*...
        ((m_max/m_star)^(n_B+1)-(m_max/m_star)^n_B)));
end
omega1 = prop.omega./...
    ((r_hat^2-prop.omega_hat)/(r_hat^2-1)+prop.r1^2*(prop.omega_hat-1)/(r_hat^2-1)/rc^2);


%-- Azimuth velocity distribution ----------------------------------------%
alpha = omega1.*((r_hat).^2-prop.omega_hat)./((r_hat).^2-1);
beta = omega1.*prop.r1.^2.*(prop.omega_hat-1)./((r_hat)^2-1);
v_theta = alpha.*r_vec+beta./r_vec;


%-- Calculate voltage and electrostatic forces ---------------------------%
V0 = m_star*prop.omega^2*rc^2*log(prop.r2/prop.r1)/e; % voltage corresponding to m_star, use single charge (q = e)
V = m_star.*log(1/r_hat)./e.*(alpha.*rc+beta./rc).^2; % voltage corresponding to m_star, use single charge (q = e)
C0 = V.*q./log(1/r_hat);
F_e = C0./r_vec; % electrostatic force


%-- Evaluate particle mobility -------------------------------------------%
if isempty(d)
    B = kernel.mp2zp(m,1,prop.T,prop.p);
else
    B = kernel.dm2zp(d,1,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B).*z; % diffusion coefficient for each mobility


%-- Axial velocity distribution ------------------------------------------%
A = pi*(prop.r2^2-prop.r1^2); % cross sectional area of APM
v_bar = prop.Q/A; % average flow velocity
v_z = 3/2*v_bar.*(1-((r_vec-rc)./(del)).^2); % axial velocity distribution
% v_z = v_bar.*ones(size(r_vec));


%-- Estimate equilibrium radius ------------------------------------------%
%-- Estimated as a vector of m and d
if round((sqrt(C0./m_star)-sqrt(C0./m_star-4*alpha*beta))/(2*alpha),15)==rc
    rs = real((sqrt(C0./m)-sqrt(C0./m-4*alpha*beta))./(2*alpha)); % equiblirium radius for a given mass
else
    rs = real((sqrt(C0./m)+sqrt(C0./m-4*alpha*beta))./(2*alpha)); % equiblirium radius for a given mass
end
lam = 2.*tau.*(alpha^2-beta^2./(rs.^4)).*prop.L./v_bar;


cond0 = or(and(rs+(prop.r2-rs).*exp(-lam)<(prop.r1-0.75*del),...
    rs+(prop.r1-rs).*exp(-lam)<(prop.r1-0.75*del)),...
    and(rs+(prop.r2-rs).*exp(-lam)>(prop.r2+6*del),...
    rs+(prop.r1-rs).*exp(-lam)>(prop.r2+6*del)));
        % NOTE: conditions limits consideration of those regions where particles
        % do not escape
ind = 1:length(m);
ind = ind(~cond0);

%-- Loop over mass, but not m_star ---------------------------------------%
Lambda = zeros(length(m),1);
for ii=ind % loop over mass, not m_star
    
    F_c = m(ii).*v_theta.^2./r_vec; % centriputal force
    c_r = tau(ii)/m(ii).*(F_c-F_e); % particle velocity across streamlines
    drcr_dr = tau(ii).*(2*alpha^2.*r_vec-2*beta^2./(r_vec.^3));
    % dcr_dr = tau/m(ii).*(...
    %     m(ii).*(alpha^2-2*alpha*beta./(r_vec.^2)-3*beta^2./(r_vec.^4))...
    %     +(q*V/log(prop.r2/prop.r1))./(r_vec.^2));

    ind_mid = 2:(length(r_vec)-1);

    zet = v_z./dz;
    gam = D(ii)/(2*dr^2);
    kap = D(ii)./(4.*r_vec.*dr);
    eta = 1./(2.*r_vec).*drcr_dr;
    % eta = 1/2.*dcr_dr;
    phi = c_r./(4*dr);

    RHS1 = @(n) (zet(1)-2*gam-eta(1)).*n(1)+...
        (gam+kap(1)-phi(1)).*n(2);
    RHS2 = @(n) (zet(ind_mid)-2*gam-eta(ind_mid)).*n(2:(end-1))+...
        (gam-kap(ind_mid)+phi(ind_mid)).*n(1:(end-2))+...
        (gam+kap(ind_mid)-phi(ind_mid)).*n(3:end);
    RHS3 = @(n) (zet(end)-2*gam-eta(end)).*n(end)+...
        (gam-kap(end)+phi(end)).*n(end-1);
    RHS = @(n) [RHS1(n),RHS2(n),RHS3(n)];

    b1 = (zet(1)+2*gam+eta(1)); % center
    c1 = (-gam-kap(1)+phi(1)); % +1

    a2 = (-gam+kap(ind_mid)-phi(ind_mid)); % -1
    b2 = (zet(ind_mid)+2*gam+eta(ind_mid)); % center
    c2 = (-gam-kap(ind_mid)+phi(ind_mid)); % +1

    a3 = (-gam+kap(end)-phi(end)); % -1
    b3 = (zet(end)+2*gam+eta(end)); % center

    a = [a2,a3];
    b = [b1,b2,b3];
    c = [c1,c2];

    n_vec0 = ones(size(r_vec));
    n_vec = n_vec0;

    % n_mat = zeros(nz,length(r_vec));
    % n_mat(1,:) = n_vec0;
    for jj = 2:nz
        n_vec = max(kernel.tridiag([0,a],b,c,RHS(n_vec)),0);
        % n_vec = max(olfert.tridiag1(a,b,c,RHS(n_vec)),0); % DON'T UNCOMMENT
        % n_mat(jj,:) = n_vec;

        Lambda_j = sum(n_vec.*v_z)/sum(n_vec0.*v_z); % updated to include v_z
        if Lambda_j < 0.0005
            Lambda_j = 0;
        end
    end
    % n_mat = max(n_mat,0);
    
    % load('L1_cmap.mat');
    % figure(1);
    % imagesc(r_vec,z_vec,n_mat);
    % colormap(cm);

    % load('L1_cmap.mat');
    % cm = 1-bsxfun(@times,cm,1-[19,19,99]./256);
    % figure(1);
    % contourf(r_vec,z_vec,n_mat,15,'edgecolor','none');
    % colormap(cm);

    Lambda(ii) = Lambda_j; % updated to include v_z
end

% figure(3);
% hold on;
% plot(m,Lambda);
% hold off;

end
