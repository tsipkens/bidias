function [B]=dp2mobility_v2(dp,P,T)

%This function converts mobility to particle size.
% dp = the particle mobility diameter (m)
% P = the pressure in atm
% T = the temperature in K
for i=1:length(dp)
    
    mu=1.81e-5*(T/293)^(0.74); %example after eq.2.29 in Hinds 
    B(i)=olfert.Cc_v3(dp(i),P,T)/3/pi/mu./dp(i);
    
end