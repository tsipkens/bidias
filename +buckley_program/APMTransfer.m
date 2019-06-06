function [tfer] = APMTransfer(z, m, invZpStar, mpStar)
% inputs to this function are the integer charge state (z) an array of
% masses (m), the mobility diameter selected by the upstream DMA (dpStar),
% and the mass corresponding to the operating conditions of the APM (mpStar)
% this transfer function calculated based off Ehara et al (1996)

clear tfer

e = 1.6022E-19; % electron charge [C]

%% APM Settings (USER INPUT)
% Radii (R2 & R1), Length (L), Rotational Speed (RPM), Aerosol Flowrate
% (Qae)

R2 = 0.025; % [m]
R1 = 0.024; % [m]
L = 0.1;    % [m]
RPM = 13350; % rotational speed in revolutions/minute
w = RPM*2*pi/60; % rotational speed in radians/second
Qae = 1.02E-3/60; % aerosol flowrate [m^3/s]

%% Transfer Function Calculation:

v_ave = Qae/(pi*(R2^2 - R1^2)); % average axial flow through APM

r_c = (R1 + R2)/2;

delAPM = (R2 - R1)/2;

V = mpStar.*r_c^2*w^2*log(R2/R1)./e;

r = (V./((m/(z*e))*w^2*log(R2/R1))).^0.5;

rho = (r-r_c)/delAPM;

tau = m/(z*e*invZpStar);

kappa = 2*tau*w^2*L/v_ave;


tfer1X = [rho > -1].*[rho < 1];
tfer2X = [rho > 1].*[rho < coth(kappa/2)];
tfer3X = [rho < -1].*[rho > -coth(kappa/2)];
tfer4X = [rho < -coth(kappa/2)];
tfer5X = [rho > coth(kappa/2)];

tfer1 = exp(-kappa);
tfer2 = ((1 - rho) + (1 + rho).*exp(-kappa))./2;
tfer3 = ((1 + rho) + (1 - rho).*exp(-kappa))./2;
% tfer4 = zeros(length(m),1);
% tfer5 = zeros(length(m),1);

tfer = tfer1X.*tfer1 + tfer2X.*tfer2 + tfer3X.*tfer3;% + tfer4X.*tfer4 + tfer5X.*tfer5;