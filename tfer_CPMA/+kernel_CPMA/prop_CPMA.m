function [prop] = prop_CPMA(opt)
%PROP_CPMA Summary of this function goes here
%   Detailed explanation goes here

if ~exist('opt','var') % if properties set is not specified
    opt = 'Buckley';
elseif isempty(opt)
    opt = 'Buckley';
end

if strcmp(opt,'Olfert')
    %-- CPMA parameters from Olfert lab --------------%
    prop.r1 = 0.06; % inner electrode radius [m]
    prop.r2 = 0.061; % outer electrode radius [m]
    prop.L = 0.2; % length of chamber [m]
    prop.p = 1; % pressure [atm]
    prop.T = 293; % system temperature [K]
    prop.Q = 3/1000/60;%0.3/1000/60;%1.5/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
    prop.omega_hat = 32/33; % ratio of angular speeds

elseif strcmp(opt,'Buckley')
    %-- CPMA/APM parameters from Buckley et al. -------------%
    prop.r2 = 0.025; % outer electrode radius [m]
    prop.r1 = 0.024; % inner electrode radius [m]
    prop.L = 0.1;    % length of APM [m]
    RPM = 13350; % rotational speed [rpm]
    prop.omega = RPM*2*pi/60; % rotational speed [rad/s]
    prop.omega_hat = 1; % APM, so rotational speed is the same
    prop.Q = 1.02e-3/60; % aerosol flowrate [m^3/s]
    prop.T = 298; % system temperature [K]
    prop.p = 1; % system pressure [atm]

end


%-- Parameters related to CPMA geometry ----------------------------------%
prop.rc = (prop.r2+prop.r1)/2;
prop.r_hat = prop.r1/prop.r2;
prop.del = (prop.r2-prop.r1)/2; % half gap width

prop.A = pi*(prop.r2^2-prop.r1^2); % cross sectional area of APM
prop.v_bar = prop.Q/prop.A; % average flow velocity


%-- For diffusion --------------------------------------------------------%
kB = 1.3806488e-23; % Boltzmann's constant
% prop.D = kB*prop.T*2.1761e+12; % diffusion coefficient
prop.D = @(B) kB.*prop.T.*B; % diffusion coefficient

end

