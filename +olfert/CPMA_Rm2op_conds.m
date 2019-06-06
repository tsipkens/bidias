function [omega,V,m_max] = CPMA_Rm2op_conds(Rm,CPMA_dims,mc,Q,q,P,T,mass_mob_pref,mass_mob_exp)

m_max = mc*(1/Rm+1);

%%%%%%calculate n_B %%%%%%%%
m_high = mc*1.001;
m_low  = mc*.999;
Bc = olfert.mass2mob(mc,mass_mob_pref,mass_mob_exp,P,T);
B_high = olfert.mass2mob(m_high,mass_mob_pref,mass_mob_exp,P,T);
B_low = olfert.mass2mob(m_low,mass_mob_pref,mass_mob_exp,P,T);
n_B = log10(B_high/B_low)/log10(m_high/m_low);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = sqrt(Q/(mc*Bc*2*pi*CPMA_dims.rc^2*CPMA_dims.L*((m_max/mc)^(n_B+1)-(m_max/mc)^n_B)));
V = mc*omega^2*CPMA_dims.rc^2*log(CPMA_dims.r2/CPMA_dims.r1)/q;





