
% PROP_DMA  Generates the prop struct used to summarize DMA parameters.
% Author:   Timothy Sipkens, 2019-07-17
%=========================================================================%

function [prop] = prop_DMA(opt_solver)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('opt_solver','var') % if order not specified
    opt_solver = 'fullydeveloped';
elseif isempty(opt_solver)
    opt_solver = 'fullydeveloped';
end


%-- DMA parameters from Olfert lab (TSI 3081) ----------------------------%
%-- +FlareNet 2018
prop.Q_s = 0.3/60/1000; % Sample flow [m^3/s]
prop.Q_a = 0.3/60/1000; % Aerosol flow [m^3/s]
prop.Q_c = 3/60/1000; % Sheath flow [m^3/s]
prop.Q_m = 3/60/1000; % Exhaust flow [m^3/s]
prop.L = 0.44369; %length, m
prop.R1 = 0.00937; %inner radius, m
prop.R2 = 0.01961; %outer radius, m
prop.T = 293; %temperature, K
prop.p = 1;

%{
%-- DMA parameters from Buckley et al. (TSI 3081) ------------------------%
prop.Q_c = 4.89E-3/60; % sheath flow [m^3/s]
prop.Q_a = 1.02E-3/60; % aerosol flow [m^3/s]
prop.Q_s = Q_a; % sample flow [m^3/s]
prop.Q_m = Q_c; % exhaust flow [m^3/s]
prop.L = 0.4444; % length of chamber [m]
prop.R2 = 0.01961; % outer electrode radius [m]
prop.R1 = 0.00937; % inner electrode radius [m]
prop.T = 294; % temperature [K]
prop.p = 1;
%}


prop.bet = (prop.Q_s+prop.Q_a)/(prop.Q_c+prop.Q_m);
prop.del = (prop.Q_s-prop.Q_a)/(prop.Q_s+prop.Q_a);

gam = (prop.R1/prop.R2)^2;
kap = prop.L/prop.R2; % Stolzenburg Manuscript Eqn 9
% kap = L*R2/(R2^2-R1^2); % Buckey et al.

%-- Calculate G_DMA ------------------------------------------------------%
switch opt_solver
    case 'Buckley' % evaluation from Buckley et al.
        I_gamma = (0.25*(1-gam^2)*(1-gam)^2+(5/18)*(1-gam^3)*(1-gam)*log(gam) + ...
            (1/12)*(1-gam^4)*(log(gam))^2)/((1-gam)*(-0.5*(1+gam)*log(gam)-(1-gam))^2);
        prop.G_DMA = 4*(1+prop.bet)^2/(1-gam)*(I_gamma+(2*(1+prop.bet)*kap)^(-2));

    case 'plug' % Plug flow, Stolzenburg, 2018
        omega_a = 1-prop.bet*(1-prop.del)/(2*(1+prop.bet))*(1-gam);
        omega_s = gam+prop.bet*(1+prop.del)/(2*(1+prop.bet))*(1-gam);

        I_gamma = @(omega) (1-omega^2)/(2*(1-gam));
        G_o_corr = (4*(1+prop.bet)^2/(1-gam)*(I_gamma(omega_s)-I_gamma(omega_a))+...
            (omega_a-omega_s)/kap^2);
        prop.G_DMA = G_o_corr*(-log(gam)/2);

    case {'Olfert','fullydeveloped'} % Olfert laboratory; Fully developed flow, Stolzenburg, 2018
        fun = @(omega) ((1-omega)^2*log(gam)/2/(1-gam)+omega*log(omega)+(1-omega))/((1+gam)*log(gam)/2+(1-gam))-prop.bet*(1-prop.del)/2/(1+prop.bet);
        guess = 0.9;
        omega_a = fzero(fun,guess);

        fun = @(omega) (1-((1-omega)^2*log(gam)/2/(1-gam)+omega*log(omega)+(1-omega))/((1+gam)*log(gam)/2+(1-gam)))-prop.bet*(1+prop.del)/2/(1+prop.bet);
        guess = 0.3;
        omega_s = fzero(fun,guess);

        A = (-1/2*(1+gam)*log(gam)-(1-gam))^(-1);
        I_gamma = @(omega) A^2/(1-gam)*...
            (-omega^2*((1-gam)*log(omega)-(1-omega)*log(gam))^2/2+...
            (omega^2*(1-gam)/2+omega^3*log(gam)/3)*((1-gam)*log(omega)-...
            (1-omega)*log(gam))+(1-omega^2)*(1-gam)^2/4+...
            5*(1-omega^3)*(1-gam)*log(gam)/18+(1-omega^4)*(log(gam))^2/12);

        G_o_corr = (4*(1+prop.bet)^2)*(I_gamma(omega_s)-I_gamma(omega_a))/(1-gam)...
            +(omega_a-omega_s)/kap^2; % Stolzenburg Manuscript Equation 8 and 12
        prop.G_DMA = G_o_corr*log(prop.R2/prop.R1); % Stolzenburg Manuscript Equation 23

    otherwise
        disp('Invalid solver specified.');
        return;

end


end

