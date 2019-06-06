function rdot = CPMA_ode_fcn(z,r_ode,v_flow,r2,r1,omega1,omega_hat,B,m,q,V)

vz = v_flow * 3/2 * (1-((2*r_ode-r2-r1)/(r2-r1))^2);
r_hat = r1/r2;
vt = omega1 * (r_hat^2-omega_hat)/(r_hat^2-1) * r_ode + omega1 * r1^2 * (omega_hat-1)/(r_hat^2-1) /r_ode;

rdot = [B/vz * (m*vt^2/r_ode - q*V/r_ode/log(r2/r1))];