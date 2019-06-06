%% Cunningham Slip Correction Factor

function ans = Cc(dp)

A1 = 1.257;
A2 = 0.4;
A3 = 0.55;

mfp = 66.5E-9;
Kn = 2*mfp./dp;
ans = 1 + Kn.*(A1 + A2*exp(-2*A3./Kn));