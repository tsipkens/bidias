
function [Lambda,prop,n] = tfer_FD(m_star,m,d,z,prop,varargin)
% TFER_FD Evaluates the transfer function of a PMA using finite differences.
% Author: Timothy Sipkens, 2018-12-27
% Based on: Buckley et al., J. Aerosol Sci. (2008) and Olfert group
% 
%-------------------------------------------------------------------------%
% Inputs:
%   m_star      Setpoint particle mass
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length)
%   varargin    Name-value pairs for setpoint    (Optional, default Rm = 3)
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Outputs:
%   Lambda      Transfer function
%   prop        Device properties (e.g. classifier length)
%   n           Struct containing information about the particle
%               distribution at different axial position (used for plotting)
%-------------------------------------------------------------------------%

if ~exist('prop','var') % if particle mass analyzer properties not specified
    prop = tfer_PMA.prop_CPMA; % will use default in this function
elseif isempty(prop)
    prop = tfer_PMA.prop_CPMA;
end

tfer_PMA.get_setpoint; % get setpoint


%-- Discretize the space -------------------------------------------------%
nr = 200;
dr = (prop.r2-prop.r1)/nr;
r_vec = (prop.r1+dr/2):dr:(prop.r2-dr/2);

nz = 200;
dz = prop.L/(nz-1);
if nargout==3; n.z_vec=0:dz:prop.L; n.r_vec=r_vec; end % vector of z positions, used for plotting only


%-- Parameters related to CPMA geometry ----------------------------------%
del = (prop.r2-prop.r1)/2;
rc = (prop.r2+prop.r1)/2;


%-- Evaluate relevant radial quantities ----------------------------------%
v_theta = sp.alpha.*r_vec+sp.beta./r_vec; % azimuthal velocity distribution
F_e = C0./r_vec; % electrostatic force
v_z = 3/2*prop.v_bar.*(1-((r_vec-rc)./(del)).^2); % axial velocity distribution (parabolic)
% v_z = v_bar.*ones(size(r_vec)); % axial velocity distribution (plug)


%-- Speed computation using resolution to limit computation --------------%
ind = 1:length(m);
if isfield(sp,'Rm') % if resolution is specified, use to reduce necessary computation
    cond0 = or(m>(m_star+2.*sp.m_max),...
        m<(m_star-2.*sp.m_max));
            % NOTE: conditions limits consideration of those regions where particles
            % do not escape, speeding computation.
    ind = ind(~cond0);
end


%-- Loop over mass, but not m_star ---------------------------------------%
if nargout==3; kk = 0; end % if outputting particle distribution, initialize kk
Lambda = zeros(1,length(m));% initialize the transfer function variable
for ii=ind % loop over mass (not m_star)
    
    F_c = m(ii).*v_theta.^2./r_vec; % centriputal force
    c_r = tau(ii)/m(ii).*(F_c-F_e); % particle velocity across streamlines
    drcr_dr = tau(ii).*(2*sp.alpha^2.*r_vec-2*sp.beta^2./(r_vec.^3));

    ind_mid = 2:(length(r_vec)-1);

    %-- Get coefficients ------------------%
    zet = v_z./dz;
    gam = D(ii)/(2*dr^2);
    kap = D(ii)./(4.*r_vec.*dr);
    eta = 1./(2.*r_vec).*drcr_dr;
    phi = c_r./(4*dr);
    
    %-- Righ-hand side of eq. to be solved ---------------%
    RHS1 = @(n) (zet(1)-2*gam-eta(1)).*n(1)+...
        (gam+kap(1)-phi(1)).*n(2);
    RHS2 = @(n) (zet(ind_mid)-2*gam-eta(ind_mid)).*n(2:(end-1))+...
        (gam-kap(ind_mid)+phi(ind_mid)).*n(1:(end-2))+...
        (gam+kap(ind_mid)-phi(ind_mid)).*n(3:end);
    RHS3 = @(n) (zet(end)-2*gam-eta(end)).*n(end)+...
        (gam-kap(end)+phi(end)).*n(end-1);
    RHS = @(n) [RHS1(n),RHS2(n),RHS3(n)];
    
    %-- Form A matrix ---------------------------------------%
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
    
    %-- Initialize variables -------------------------%
    n_vec = ones(size(r_vec)); % number concentration at current axial position
    n_vec0 = n_vec; % initial number concentration (used for tfer_PMA. func. eval.)
    if nargout==3 % initilize variables used for visualizing number concentrations
        n_mat = zeros(nz,length(r_vec));
        n_mat(1,:) = n_vec;
    end
    
    %-- Primary loop for finite difference ---------------------------%
    for jj = 2:nz
        n_vec = max(tfer_PMA.tridiag([0,a],b,c,RHS(n_vec)),0);
            % solve system using Thoman algorithm
            
        if nargout==3; n_mat(jj,:) = n_vec; end
            % record particle distribution at this axial position
    end
    
    if nargout==3 % for visualizing number concentrations
        kk = kk+1;
        n.n_mat{kk} = max(n_mat,0);
    end
    
    Lambda(ii) = sum(n_vec.*v_z)/sum(n_vec0.*v_z);
        % evaluate transfer fucntion
        
    if Lambda(ii) < 0.0005; Lambda(ii) = 0; end
        % truncate small tfer_PMA. func. values
    
end

end
