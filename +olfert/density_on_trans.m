close all;
clear;

%mass_mob_pref=524;
mass_mob_pref=[1130 524 524/2];
mass_mob_exp = [3 3 3];
%omega = 702; 
%V       = 186;
omega_hat = .95;
q       =   1 * 1.6 *10^-19;    %particle charge (C) assume one charge per particle (boltzmann dist.)
Q       =   0.3/1000/60;        %volume flow rate (m^3/s) ~1 lpm
T       =   293;                %temperature (K) 
P       =   1;                  %pressure in (atm)

[CPMA_dims] = CPMA_Camb_dimstest(omega_hat);   %load Cambustion CPMA dimensions

Rm =5
mc = 0.1e-18;

%[Rm,mc,m_max]=CPMA_Camb_Res(mass_mob_pref,mass_mob_exp,omega,q,V,Q,T,P,CPMA_dims)

for j=1:length(mass_mob_exp)
[omega,V,m_max] = CPMA_Rm2op_conds(Rm,CPMA_dims,mc,Q,q,P,T,mass_mob_pref(2),mass_mob_exp(2));

omega
V

x=2*mc-m_max;
if x< 0
    m_min = 0;
else
    m_min=2*mc-m_max;
end

m(:,j) = linspace(m_min/1.1,m_max*1.1,200);
B(:,j) = mass2mob(m(:,j),mass_mob_pref(j),mass_mob_exp(j),P,T);


for i=1:length(m)

[trans(i,j)]=CPMA_nodiff_model(m(i,j),B(i,j),q,Q,omega,V,CPMA_dims);
%[trans_D(i)]=CPMA_diff_model(m(i),B(i),q,Q,omega,V,T,CPMA_dims);
%[trans(i)]=CPMA_diff_model_fixedgrid(m(i,j),B(i,j),q,Q,omega,V,T,CPMA_dims);
%[trans(i,j)]=triangular_transfer_fcn(1/Rm,m(i)/mc);

end

end
 figure
 plot(m(:,1)*1e18,trans(:,1),m(:,2)*1e18,trans(:,2),m(:,3)*1e18,trans(:,3));
 legend('NaCl','H2O','"soot"');
 xlabel('mass (fg)');
ylabel('transfer function');
 
nacl_error=trapz(m(:,1),trans(:,1))/trapz(m(:,2),trans(:,2))
 
soot_error=trapz(m(:,3),trans(:,3))/trapz(m(:,2),trans(:,2))
 

r_hat_sq = CPMA_dims.r_hat^2
omega_hat=omega_hat

if CPMA_dims.r_hat^2 > omega_hat
    disp('Taylor limit exceeded');
else
    disp('Taylor limit not exceeded');
end

 %figure
 %plot(m,trans_Df-trans_D);
