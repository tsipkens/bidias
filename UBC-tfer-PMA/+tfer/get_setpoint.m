
% GET_SETPOINT Evaluates setpoint parameters, including C0, alpha, and beta.
% Author: Timothy A. Sipkens, 2019-05-01
% Note: As a script, this code uses variables currently in the workspace. 
%  This script is also used to parse some of the inputs to the various 
%  transfer functions, including the existence of the integer charge state 
%  and particle mobility. 
%
%-------------------------------------------------------------------------%
% Requied variables:
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
% Smaple outputs:
%   C0          Summary parameter for the electrostatic force
%   tau         Product of mechanical mobility and particle mass
%   sp          Struct containing mutliple setpoint parameters (V, alpha, etc.)
%-------------------------------------------------------------------------%


%-- Set up mobility calculations and parse inputs ------------------------%
if ~exist('z','var') % if integer charge is not specified, use z = 1
    z = 1;
end
e = 1.60218e-19; % electron charge [C]
q = z.*e; % particle charge

if ~exist('d','var') % evaluate mechanical mobility
    warning('Invoking mass-mobility relation to determine Zp.');
    B = tfer.mp2zp(m,z,prop.T,prop.p);
else
    B = tfer.dm2zp(d,z,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B).*z; % diffusion as a function of mechanical mobiltiy and charge state


%-- Parse inputs for setpoint --------------------------------------------%
if exist('varargin','var') % parse input name-value pair, if it exists
    if length(varargin)==2
        sp.(varargin{1}) = varargin{2};
    else
        sp.Rm = 3; % by default use resolution, with value of 3
    end
else
    sp.Rm = 3; % by default use resolution, with value of 3
end


%-- Proceed depnding on which setpoint parameter is specified ------------%
if isfield(sp,'omega1') % if angular speed of inner electrode is specified
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    sp.V = m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
    
    
elseif isfield(sp,'V') % if voltage is specified
    
    v_theta_rc = sqrt(sp.V.*e./(m_star.*log(1/prop.r_hat)));
    A = prop.rc.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1)+...
        1./prop.rc.*(prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1));
    sp.omega1 = v_theta_rc./A;
    
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    
elseif isfield(sp,'Rm') % if resolution is specified
    
    %-- Use definition of Rm to derive angular speed at centerline -------%
    %-- See Reavell et al. (2011) for resolution definition --%
    n_B = -0.6436;
    B_star = tfer.mp2zp(m_star,1,prop.T,prop.p); % involves invoking mass-mobility relation
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
    sp.V = m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
    
    
else
    error('Invalid setpoint parameter specified.');
    
end


C0 = sp.V.*q./log(1/prop.r_hat); % calcualte recurring C0 parameter


