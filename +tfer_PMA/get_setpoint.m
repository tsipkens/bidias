
% GET_SETPOINT  Script to evaluate setpoint parameters including C0, alpha, and beta.
% Author:       Timothy A. Sipkens, 2019-05-01
%=========================================================================%
%
%-------------------------------------------------------------------------%
% Required variables:
%   m_star      Mass corresponding to the measurement set point of the APM
%   d           Particle mobility diameter, can be vector [nm]
%   m           Particle mass, can be vector of same length as d
%   z           Integer charge state, scalar
%   prop        Properties of particle mass analyzer
%   varargin    Name-value pairs for setpoint    (Optional, default Rm = 3)
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Sample outputs:
%   C0          Summary parameter for the electrostatic force
%   tau         Product of mechanical mobility and particle mass
%   sp          Struct containing mutliple setpoint parameters (V, alpha, etc.)
% 
% Notes:    As a script, this code uses variables currently in the 
%           workspace. This script is also used to parse some of the inputs 
%           to the various transfer functions, including the existence of 
%           the integer charge state and particle mobility. 
%-------------------------------------------------------------------------%


%-- Parse inputs ---------------------------------------------------------%
if ~exist('z','var'); z = []; end
if isempty(z); z = 1; end % if integer charge is not specified, use z = 1

if ~exist('m_star','var'); m_star = []; end


%-- Parse inputs for setpoint --------------------------------------------%
if exist('varargin','var') % parse input name-value pair, if it exists
    if length(varargin)==2
        sp.(varargin{1}) = varargin{2};
    elseif length(varargin)==4
        sp.(varargin{1}) = varargin{2};
        sp.(varargin{3}) = varargin{4};
    else
        sp.Rm = 3; % by default use resolution, with value of 3
    end
else
    sp.Rm = 3; % by default use resolution, with value of 3
end


%-- Set up mobility calculations -----------------------------------------%
e = 1.60218e-19; % electron charge [C]
q = z.*e; % particle charge

if ~exist('d','var') % evaluate mechanical mobility
    warning('Invoking mass-mobility relation to determine Zp.');
    B = tfer_PMA.mp2zp(m,z,prop.T,prop.p);
else
    B = tfer_PMA.dm2zp(d,z,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B).*z; % diffusion as a function of mechanical mobiltiy and charge state


%=========================================================================%
%-- Proceed depending on which setpoint parameter is specified -----------%
if isempty(m_star) % case if m_star is not specified (use voltage and speed)
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega.*prop.rc.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    m_star = sp.V./(log(1/prop.r_hat)./e.*...
        (sp.alpha.*prop.rc+sp.beta./prop.rc).^2);
        % q = e, z = 1 for setpoint
    
    sp.omega1 = sp.alpha+sp.beta./(prop.r1.^2);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    
elseif isfield(sp,'omega1') % if angular speed of inner electrode is specified
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    sp.V = m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
        % q = e, z = 1 for setpoint
    
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    sp.omega = sp.alpha+sp.beta./(prop.rc.^2);
        
elseif isfield(sp,'omega') % if angular speed at gap center is specified
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega.*prop.rc.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    sp.V = m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
        % q = e, z = 1 for setpoint
    
    sp.omega1 = sp.alpha+sp.beta./(prop.r1.^2);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    
elseif isfield(sp,'V') % if voltage is specified
    
    v_theta_rc = sqrt(sp.V.*e./(m_star.*log(1/prop.r_hat)));
        % q = e, z = 1 for setpoint
    A = prop.rc.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1)+...
        1./prop.rc.*(prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1));
    sp.omega1 = v_theta_rc./A;
    
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    sp.omega = sp.alpha+sp.beta./(prop.rc.^2);
    
elseif isfield(sp,'Rm') % if resolution is specified
    
    %-- Use definition of Rm to derive angular speed at centerline -------%
    %-- See Reavell et al. (2011) for resolution definition --%
    n_B = -0.6436;
    B_star = tfer_PMA.mp2zp(m_star,1,prop.T,prop.p);
        % involves invoking mass-mobility relation
        % z = 1 for the setpoint
        
    sp.m_max = m_star*(1/sp.Rm+1);
    omega = sqrt(prop.Q/(m_star*B_star*2*pi*prop.rc^2*prop.L*...
        ((sp.m_max/m_star)^(n_B+1)-(sp.m_max/m_star)^n_B)));
    
    %-- Evaluate angular speed of inner electrode ------------------------%
    sp.omega1 = omega./...
        ((prop.r_hat^2-prop.omega_hat)/(prop.r_hat^2-1)+...
        prop.r1^2*(prop.omega_hat-1)/(prop.r_hat^2-1)/prop.rc^2);
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    sp.omega = sp.alpha+sp.beta./(prop.rc.^2);
    sp.V = m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
    
else
    error('Invalid setpoint parameter specified.');
    
end
sp.m_star = m_star;
    % copy center mass to sp
    % center mass is for a singly charged particle

B_star = tfer_PMA.mp2zp(m_star,1,prop.T,prop.p);
C0 = sp.V.*q./log(1/prop.r_hat); % calcualte recurring C0 parameter


