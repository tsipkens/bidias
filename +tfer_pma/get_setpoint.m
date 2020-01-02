
% GET_SETPOINT Generate setpoint parameter structure from available parameters. 
% Author:  Timothy Sipkens, 2019-11-21
%-------------------------------------------------------------------------%
% Required variables:
%   prop        Properties of particle mass analyzer
%   varargin    Name-value pairs for setpoint    (Optional, default Rm = 3)
%                   ('Rm',double) - Resolution
%                   ('omega1',double) - Angular speed of inner electrode
%                   ('V',double) - Setpoint voltage
%
% Sample outputs:
%   sp          Struct containing mutliple setpoint parameters (V, alpha, etc.)
%   m_star      Setpoint mass, assuming a singly charged particle
% 
% Notes:    As a script, this code uses variables currently in the 
%           workspace. This script is also used to parse some of the inputs 
%           to the various transfer functions, including the existence of 
%           the integer charge state and particle mobility. 
%=========================================================================%

function [sp,m_star] = get_setpoint(prop,varargin)


%-- Initialize sp structure ----------------------------------------------%
sp = struct('m_star',[],'V',[],'Rm',[],'omega',[],...
    'omega1',[],'omega2',[],'alpha',[],'beta',[],'m_max',[]);
        % default empty structure
        % sp.m_star corresponds the mass at the given setpoint
        %   for a singly-charged particle


%-- Parse inputs ---------------------------------------------------------%
if length(varargin)==2 % only a single setpoint parameter is specified, assume Rm = 3
    [sp.Rm] = 3; % by default use resolution, with value of 3
else
    for ii=1:2:length(varargin)
        sp.(varargin{ii}) = varargin{ii+1};
    end
end
%-------------------------------------------------------------------------%


e = 1.60218e-19; % electron charge [C]


%== Proceed depending on which setpoint parameters are specified =========%
%== CASE 1: m_star is not specified, use V and omega =====================%
if isempty(sp.m_star)
    % case if m_star is not specified
    % requires that voltage, 'V', and speed, 'omega' are specified
    
    %-- Evaluate angular speed of inner electrode ------------------------%
    if isempty(sp.omega1)
        sp.omega1 = sp.omega./...
            ((prop.r_hat^2-prop.omega_hat)/(prop.r_hat^2-1)+...
            prop.r1^2*(prop.omega_hat-1)/(prop.r_hat^2-1)/prop.rc^2);
    end
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    sp.m_star = sp.V./(log(1/prop.r_hat)./e.*...
        (sp.alpha.*prop.rc+sp.beta./prop.rc).^2);
        % q = e, z = 1 for setpoint
    
    sp.omega = sp.alpha+sp.beta./(prop.rc.^2);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    
    
%== CASE 2: m_star and omega1 are specified ==============================%
elseif ~isempty(sp.omega1) % if angular speed of inner electrode is specified
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    sp.V = sp.m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
        % q = e, z = 1 for setpoint
    
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    sp.omega = sp.alpha+sp.beta./(prop.rc.^2);
    
    
%== CASE 3: m_star and omega are specified ===============================%
elseif ~isempty(sp.omega) % if angular speed at gap center is specified
    
    %-- Evaluate angular speed of inner electrode ------------------------%
    sp.omega1 = sp.omega./...
        ((prop.r_hat^2-prop.omega_hat)/(prop.r_hat^2-1)+...
        prop.r1^2*(prop.omega_hat-1)/(prop.r_hat^2-1)/prop.rc^2);
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    
    sp.V = sp.m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
        % q = e, z = 1 for setpoint
    
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    
    
%== CASE 4: m_star and V are specified ===================================%
elseif ~isempty(sp.V) % if voltage is specified
    
    v_theta_rc = sqrt(sp.V.*e./(sp.m_star.*log(1/prop.r_hat)));
        % q = e, z = 1 for setpoint
    A = prop.rc.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1)+...
        1./prop.rc.*(prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1));
    sp.omega1 = v_theta_rc./A;
    
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    sp.omega = sp.alpha+sp.beta./(prop.rc.^2);
    

%== CASE 5: m_star and Rm are specified ==================================%
elseif ~isempty(sp.Rm) % if resolution is specified
    
    %-- Use definition of Rm to derive angular speed at centerline -------%
    %-- See Reavell et al. (2011) for resolution definition --%
    n_B = -0.6436;
    B_star = tfer_pma.mp2zp(sp.m_star,1,prop.T,prop.p,prop);
        % involves invoking mass-mobility relation
        % z = 1 for the setpoint
    
    sp.m_max = sp.m_star*(1/sp.Rm+1);
    sp.omega = sqrt(prop.Q/(sp.m_star*B_star*2*pi*prop.rc^2*prop.L*...
        ((sp.m_max/sp.m_star)^(n_B+1)-(sp.m_max/sp.m_star)^n_B)));
    
    %-- Evaluate angular speed of inner electrode ------------------------%
    sp.omega1 = sp.omega./...
        ((prop.r_hat^2-prop.omega_hat)/(prop.r_hat^2-1)+...
        prop.r1^2*(prop.omega_hat-1)/(prop.r_hat^2-1)/prop.rc^2);
    
    %-- Azimuth velocity distribution and voltage ------------------------%
    sp.alpha = sp.omega1.*(prop.r_hat.^2-prop.omega_hat)./(prop.r_hat.^2-1);
    sp.beta = sp.omega1.*prop.r1.^2.*(prop.omega_hat-1)./(prop.r_hat^2-1);
    sp.omega2 = sp.alpha+sp.beta./(prop.r2.^2);
    sp.V = sp.m_star.*log(1/prop.r_hat)./e.*(sp.alpha.*prop.rc+sp.beta./prop.rc).^2;
    
else % other cases are not supported, output error
    error('Invalid setpoint parameters specified.');
    
end


%-- Calculate resolution -------------------------------------------------%
if isempty(sp.Rm) % if resolution is not specified
    n_B = -0.6436;
    B_star = tfer_pma.mp2zp(sp.m_star,1,prop.T,prop.p,prop);
    t0 = prop.Q/(sp.m_star*B_star*2*pi*prop.L*...
        sp.omega^2*prop.rc^2);
    m_rat = @(Rm) 1/Rm+1;
    fun = @(Rm) (m_rat(Rm))^(n_B+1)-(m_rat(Rm))^n_B;
    sp.Rm = fminsearch(@(Rm) (t0-fun(Rm))^2,5);
    sp.m_max = sp.m_star*(1/sp.Rm+1);
end


m_star = sp.m_star; % output m_star independently

end
