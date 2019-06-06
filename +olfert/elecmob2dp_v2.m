function [dp]=elecmob2dp(Zp,n,P,T)

%This function converts electrical mobility to particle size.
% Zp = the particle electrical mobility
% n = the charge on the particle
% P = the pressure in atm
% T = the temperature in K

e = 1.602e-19;      %electron charge
mu=1.81e-5*(T/293)^(0.74); %example after eq.2.29 in Hinds 

for i=1:n
    for j=1:length(Zp)
        f=@(dp)i*e/3/pi/mu/Zp(j)-dp/olfert.Cc_v3(dp,P,T); % Orig: CC_v2
        
        %guess a solution from a fit i did
        guess=(7.571e-12)*Zp(j)^(-5.569e-1);
        
        dp(j,i) = fzero(f,guess);
    end
end