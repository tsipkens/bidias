function [Rm,mc,m_max] = CPMA_Camb_Res(mass_mob_pref,mass_mob_exp,omega,q,V,Q,T,P,CPMA_dims)



mc = q*V/omega^2/CPMA_dims.rc^2/log(CPMA_dims.r2/CPMA_dims.r1);
Bc = mass2mob(mc,mass_mob_pref,mass_mob_exp,P,T);



f=@(m_max)omega^2*CPMA_dims.rc*m_max-q*V/CPMA_dims.rc/log(CPMA_dims.r2/CPMA_dims.r1)-Q/2/pi/mass2mob(m_max,mass_mob_pref,mass_mob_exp,P,T)/CPMA_dims.L/CPMA_dims.rc;
guess=(Q/2/pi/Bc/CPMA_dims.L/CPMA_dims.rc+q*V/CPMA_dims.rc/log(CPMA_dims.r2/CPMA_dims.r1))/(omega^2*CPMA_dims.rc);
options = optimset('TolX',mc/1000); %need to set the tolerance to for these small masses
[m_max fval exitflag output] = fzero(f,guess,options);


Rm = 1/(m_max/mc-1);


