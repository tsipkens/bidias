function [B]=dp2mobility_v3(dp,P,T)

%This function converts mobility to particle size.
% dp = the particle mobility diameter (m)
% P = the pressure in atm
% T = the temperature in K
for i=1:length(dp)
    mu=1.81809e-5*(T/293.15)^1.5*(1/P)*(293.15+110.4)/(T+110.4); %see Rader (1990)
    B(i)=olfert.Cc_v3(dp(i),P,T)/3/pi/mu./dp(i);
end