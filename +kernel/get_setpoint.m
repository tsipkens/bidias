
%-- Set up mobility calculations -----------------------------------------%
if ~exist('z','var') % if integer charge is not specified, use z = 1
    z = 1;
end
e = 1.60218e-19; % electron charge [C]
q = z.*e;

if ~exist('d','var')
    warning('Invoking mass-mobility relation to determine Zp.');
    B = kernel.mp2zp(m,z,prop.T,prop.p);
else
    B = kernel.dm2zp(d,z,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B).*z;


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
    B_star = kernel.mp2zp(m_star,1,prop.T,prop.p); % involves invoking mass-mobility relation
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


