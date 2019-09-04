
% PROP_CPMA Generates the prop struct used to summarize CPMA parameters.
% Author:   Timothy Sipkens, 2019-06-26
%=========================================================================%

function [prop] = prop_CPMA(opt)
%-------------------------------------------------------------------------%
% Input:
%   opt         Options string specifying parameter set
%                   (Optional, default 'Olfert')
%
% Output:
%   prop        Properties struct for use in evaluating transfer function
%-------------------------------------------------------------------------%


if ~exist('opt','var') % if properties set is not specified
    opt = 'Olfert';
elseif isempty(opt)
    opt = 'Olfert';
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
    prop.omega_hat = 1; % APM, so rotational speed is the same
    prop.Q = 1.02e-3/60; % aerosol flowrate [m^3/s]
    prop.T = 298; % system temperature [K]
    prop.p = 1; % system pressure [atm]

elseif strcmp(opt,'Ehara')
    %-- APM parameters from Ehara et al. -------------%
    prop.r2 = 0.103; % outer electrode radius [m]
    prop.r1 = 0.1; % inner electrode radius [m]
    prop.L = 0.2;    % length of APM [m]
    prop.omega_hat = 1; % APM, so rotational speed is the same
    prop.Q = 0.5/1000/60; % aerosol flowrate [m^3/s], assumed
    prop.T = 298; % system temperature [K]
    prop.p = 1; % system pressure [atm]

elseif strcmp(opt,'Olfert-Collings')
    %-- Parameters from Olfert and Collings -------------%
    %   Nearly identical to the Ehara et al. case
    prop.r2 = 0.103; % outer electrode radius [m]
    prop.r1 = 0.1; % inner electrode radius [m]
    prop.L = 0.2;
    prop.omega_hat = 0.945;
    prop.Q = 0.5/1000/60; % aerosol flowrate [m^3/s]
    prop.T = 295; % system temperature [K]
    prop.p = 1; % system pressure [atm]
    
elseif strcmp(opt,'Kuwata')
    %-- Parameters from Kuwata --------------------------%
    prop.r2 = 0.052; % outer electrode radius [m]
    prop.r1 = 0.05; % inner electrode radius [m]
    prop.L = 0.25;
    prop.omega_hat = 1;
    prop.Q = 1.67e-5; % aerosol flowrate [m^3/s]
    prop.T = 295; % system temperature [K]
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
prop.D = @(B) kB.*prop.T.*B; % diffusion coefficient


end

