
function [Lambda,prop,F] = tfer_CPMA_FD(m_star,m,d,z,prop,varargin)
% TFER_CPMA_FD Evaluates the transfer function for a CPMA using finite difference.
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

opt_pd = (nargout==3); % output solutions for plotting

kernel.get_setpoint; % get setpoint


%-- Discretize the space -------------------------------------------------%
nr = 200;
dr = (prop.r2-prop.r1)/nr;
r_vec = (prop.r1+dr/2):dr:(prop.r2-dr/2);

nz = 200;
dz = prop.L/(nz-1);
if opt_pd; F.z_vec=0:dz:prop.L; F.r_vec=r_vec; end % vector of z locations, used for plotting only


%-- Parameters related to CPMA geometry ----------------------------------%
del = (prop.r2-prop.r1)/2;
rc = (prop.r2+prop.r1)/2;


%-- Evaluate relevant radial quantities ----------------------------------%
v_theta = sp.alpha.*r_vec+sp.beta./r_vec; % azimuthal velocity distribution
F_e = C0./r_vec; % electrostatic force
v_z = 3/2*prop.v_bar.*(1-((r_vec-rc)./(del)).^2); % axial velocity distribution (parabolic)
% v_z = v_bar.*ones(size(r_vec)); % axial velocity distribution (plug)


%-- Estimate equilibrium radius ------------------------------------------%
%-- Estimated as a vector of m and d --%
if round((sqrt(C0./m_star)-sqrt(C0./m_star-4*sp.alpha*sp.beta))/(2*sp.alpha),15)==rc
    rs = real((sqrt(C0./m)-sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
else
    rs = real((sqrt(C0./m)+sqrt(C0./m-4*sp.alpha*sp.beta))./(2*sp.alpha)); % equiblirium radius for a given mass
end
lam = 2.*tau.*(sp.alpha^2-sp.beta^2./(rs.^4)).*prop.L./prop.v_bar;

ind = 1:length(m);
if isfield(sp,'Rm') % if resolution is specified, use to reduce necessary computation
    cond0 = or(m>(m_star+2.*sp.m_max),...
        m<(m_star-2.*sp.m_max));
            % NOTE: conditions limits consideration of those regions where particles
            % do not escape, speeding computation.
    ind = ind(~cond0);
end


%-- Loop over mass, but not m_star ---------------------------------------%
if opt_pd; kk = 0; end
Lambda = zeros(1,length(m));
for ii=ind % loop over mass, not m_star
    
    F_c = m(ii).*v_theta.^2./r_vec; % centriputal force
    c_r = tau(ii)/m(ii).*(F_c-F_e); % particle velocity across streamlines
    drcr_dr = tau(ii).*(2*sp.alpha^2.*r_vec-2*sp.beta^2./(r_vec.^3));

    ind_mid = 2:(length(r_vec)-1);

    zet = v_z./dz;
    gam = D(ii)/(2*dr^2);
    kap = D(ii)./(4.*r_vec.*dr);
    eta = 1./(2.*r_vec).*drcr_dr;
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
    
    if opt_pd % for visualizing number concentrations
        n_mat = zeros(nz,length(r_vec));
        n_mat(1,:) = n_vec0;
    end
    
    for jj = 2:nz
        n_vec = max(kernel.tridiag([0,a],b,c,RHS(n_vec)),0);
        if opt_pd; n_mat(jj,:) = n_vec; end

        Lambda_j = sum(n_vec.*v_z)/sum(n_vec0.*v_z); % updated to include v_z
        if Lambda_j < 0.0005
            Lambda_j = 0;
        end
    end
    
    if opt_pd % for visualizing number concentrations
        kk = kk+1;
        F.n_mat{kk} = max(n_mat,0);
    end

    % load('L1_cmap.mat');
    % cm = 1-bsxfun(@times,cm,1-[19,19,99]./256);

    Lambda(ii) = Lambda_j; % updated to include v_z
end

% figure(3);
% hold on;
% plot(m,Lambda);
% hold off;

end
