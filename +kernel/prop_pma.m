
% PROP_PMA  Generates the prop struct used to summarize CPMA parameters.
% Author:   Timothy Sipkens, 2019-06-26
%=========================================================================%

function [prop] = prop_pma(opts)
%-------------------------------------------------------------------------%
% Input:
%   opt         Options string specifying parameter set
%                   (Optional, default 'Olfert')
%
% Output:
%   prop        Properties struct for use in evaluating transfer function
%-------------------------------------------------------------------------%


if ~exist('opts','var') % if properties set is not specified
    opts = 'Olfert';
elseif isempty(opts)
    opts = 'Olfert';
end

prop.mass_mob_pref = 524;
prop.mass_mob_exp = 3;

switch opts
    
    %-- CPMA parameters from Olfert lab ----------------------------------%
    case 'Olfert'
        prop.r1 = 0.06; % inner electrode radius [m]
        prop.r2 = 0.061; % outer electrode radius [m]
        prop.L = 0.2; % length of chamber [m]
        prop.p = 1; % pressure [atm]
        prop.T = 293; % system temperature [K]
        prop.Q = 3/1000/60;%0.3/1000/60;%1.5/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
        prop.omega_hat = 32/33; % ratio of angular speeds

    %-- CPMA/APM parameters from Buckley et al. --------------------------%
    case 'Buckley'
        prop.r2 = 0.025; % outer electrode radius [m]
        prop.r1 = 0.024; % inner electrode radius [m]
        prop.L = 0.1;    % length of APM [m]
        RPM = 13350; % rotational speed [rpm]
        prop.omega = RPM*2*pi/60; % rotational speed [rad/s]
        prop.omega_hat = 1; % APM, so rotational speed is the same
        prop.Q = 1.02e-3/60; % aerosol flowrate [m^3/s]
        prop.T = 298; % system temperature [K]
        prop.p = 1; % system pressure [atm]
    
    %-- CPMA parameters from Olfert lab ----------------------------------%
    case {'FlareNet18','fn18'}
        prop.r1 = 0.06; % inner electrode radius [m]
        prop.r2 = 0.061; % outer electrode radius [m]
        prop.L = 0.2; % length of chamber [m]
        prop.p = 1; % pressure [atm]
        prop.T = 293; % system temperature [K]
        prop.Q = 0.3/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
        prop.omega_hat = 32/33; % ratio of angular speeds

    %-- APM parameters from Ehara et al. -------------%
    case 'Ehara'
        prop.r2 = 0.103; % outer electrode radius [m]
        prop.r1 = 0.1; % inner electrode radius [m]
        prop.L = 0.2;    % length of APM [m]
        prop.omega_hat = 1; % APM, so rotational speed is the same
        prop.Q = 0.5/1000/60; % aerosol flowrate [m^3/s], assumed
        prop.T = 298; % system temperature [K]
        prop.p = 1; % system pressure [atm]
    
    %-- Parameters from Olfert and Collings -------------%
    %   Nearly identical to the Ehara et al. case
    case 'Olfert-Collings'
        prop.r2 = 0.103; % outer electrode radius [m]
        prop.r1 = 0.1; % inner electrode radius [m]
        prop.L = 0.2;
        prop.omega_hat = 0.945;
        prop.Q = 0.5/1000/60; % aerosol flowrate [m^3/s]
        prop.T = 295; % system temperature [K]
        prop.p = 1; % system pressure [atm]
    
    %-- Parameters from Kuwata --------------------------%
    case 'Kuwata'
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

