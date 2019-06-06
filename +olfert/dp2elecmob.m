function [Zp]=dp2elecmob(dp,n,P,T)

%This function converts electrical mobility to particle size.
% dp = the particle mobility diameter (m)
% n = the number of charges on the particle
% P = the pressure in atm
% T = the temperature in K


%% mu
S = 110.4; % Temperature in K
T_0 = 296.15; % Reference temperature in K
vis_23 = 1.83245*10^-5; % Reference viscosity in kg/(m*s) or Pa*s or N*s/m2

% Kim et al. (2005), ISO 15900 Eqn 3
mu = vis_23 * ((T / T_0)^1.5)*((T_0 + S)/(T + S));

%% mfp
S = 110.4; % Temperature in K
mfp_0 = 6.730*10^-8; % Mean free path of gas molecules in air at reference conditions in m
T_0 = 296.15; % Reference temperature in K
P_0 = 101325; % Reference pressure in Pa (760 mmHg to Pa)

P = P * P_0;

% Kim et al. (2005) (10.6028/jres.110.005), ISO 15900 Eqn 4
mfp = mfp_0 * (T / T_0)^2 * (P_0 / P) * ((T_0 + S)/(T + S));

%% Cc
Kn = 2 .* mfp ./ dp;
% Kim et at. (2005) (910.6028/jres.110.005) or  ISO 15900 Table 1
alpha = 1.165;
beta = 0.483;
gamma = 0.997;
Cc = 1 + Kn .* (alpha + beta .* exp(-gamma ./ Kn));

q = 1.60217646*10^-19; %elementary charge

Zp=(n*q.*Cc)./(3*pi*mu.*dp);

