
omega_hat = 1;
[CPMA_dims] = olfert.CPMA_Camb_dims(omega_hat);

RPM = 13350; % rotational speed [rpm]
omega = RPM*2*pi/60; % rotational speed [rad/s]
Q = 1.02E-3/60; % aerosol flowrate [m^3/s]
T = 298; % system temperature [K]

m_star = 1.0;
B = 1;
q = 1*1.6022E-19; % electron charge [C];
V = 29.8902; % m_star.*r_c^2*omega^2*log(R2/R1)./e;

[T_old] = olfert.CPMA_diff_model_fixedgrid(m_star,B,q,Q,omega,V,T,CPMA_dims);






