%% Mobility Diameter from Electrical Mobility Calculation using Cunningham Slip Correction Factor (iteration)
function [dp] = zp2dpCc(Zp)
% Zp should be in m^2/(V*s)
n = 1; %assume singly charged
e = 1.6022E-19; %electron charge
mu = 1.79E-5; %dynamic viscosity
mair = 29.21*1.66053892E-27; % mass of nitrogen argon mixture for plasma measurements (from summary.xls shared doc)
% mair = 28.97*1.66053892E-27;
kb = 1.3806488E-23;
rho = 1.24; %density of gas
T = 298;
mfp = 65.8; % mean free path of gas mixture in nm (from summary.xls)

if (1/Zp) < 5E6
    eps = 1.36;
    coeff = (mair*pi/(8*kb*T))^0.5 * 3*e/(pi*rho*eps);
    dp0 = (coeff/Zp)^0.5 - 0.3E-9; % initial guess for FM regime
    dp_FM = dp0;
else
%     dp0 = n*e/(3*pi*mu*Zp) - 0.3E-9; %initial guess fo continuum regime at mobility diameter, dp.
    dp0 = n*e/(3*pi*mu*Zp);
end

% Cunningham Slip Correction Factor (currently guess at formula)
a1 = 1.257;
a2 = 0.4;
a3 = 0.55;
i = 1;
% Kn = 2*mfp*1E-9/(dp0 + 0.3E-9); %Knudsen number
Kn = 2*mfp*1E-9/(dp0); %Knudsen number
cFactor = 1 + Kn*(a1 + a2*exp(-2*a3/Kn)); %Slip correction factor
% dp = n*e*cFactor/(3*pi*mu*Zp) - 0.3E-9;
dp = n*e*cFactor/(3*pi*mu*Zp);
dpMAT(i) = dp0;
dpMAT(i+1) = dp;

while abs(dp0-dp)/dp0 > 0.0001 && i < 1000
    i = i + 1;
    dp0 = (dp + dp0)/2;
%     Kn = 2*mfp*1E-9/(dp0+0.3E-9); %Knudsen number
    Kn = 2*mfp*1E-9/(dp0); %Knudsen number
    cFactor = 1 + Kn*(a1 + a2*exp(-2*a3/Kn));
%     dp = n*e*cFactor/(3*pi*mu*Zp) - 0.3E-9;
     dp = n*e*cFactor/(3*pi*mu*Zp);
    dpMAT(i) = dp;
end
