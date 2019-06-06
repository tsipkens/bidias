function cun =   Cc_v3(d,P,T)

% This function gives the cunningham slip correction
% factor of particles with diamter d (m) at pressure P (atm) and T (K)

% Cunningham slip correction parameters (Jae Hyun Kim, George W. Mulholland, S R. Kukuck, D Y. Pui, 2005)
alpha   = 1.165*2;
beta    = 0.483*2;
gamma   = .997/2;

la      = 67.30e-9;        % mfp of air at 101.325kPa and 296.15 K
lap     = la*(T/296.15)^2*(1/P)*((110.4+296.15)/(T+110.4));% mfp of air at new P and T

cun     = 1+((lap)./d).*(alpha + beta*exp(-gamma*d/lap));
