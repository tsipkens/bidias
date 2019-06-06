
% close all;
clear;

mass_mob_pref=524;
%mass_mob_pref=1130;
mass_mob_exp = 3;
%omega = 702; 
%V       = 186;
omega_hat = 32/33;
q       =   1 * 1.6 *10^-19;    %particle charge (C) assume one charge per particle (boltzmann dist.)
Q       =   1.5/1000/60;        %volume flow rate (m^3/s) ~1 lpm
T       =   293;                %temperature (K) 
P       =   1;                  %pressure in (atm)

[CPMA_dims] = olfert.CPMA_Camb_dims(omega_hat);   %load Cambustion CPMA dimensions

Rm =10;
mc = 1e-20;
Bc=olfert.mass2mob(mc,mass_mob_pref,mass_mob_exp,P,T);
dp = olfert.elecmob2dp_v2(Bc*q,1,P,T)*1e9;

%[Rm,mc,m_max]=CPMA_Camb_Res(mass_mob_pref,mass_mob_exp,omega,q,V,Q,T,P,CPMA_dims)


[omega,V,m_max] = olfert.CPMA_Rm2op_conds(Rm,CPMA_dims,mc,Q,q,P,T,mass_mob_pref,mass_mob_exp);
omega_rpm = omega/2/pi*60;
x=2*mc-m_max;
if x< 0
    m_min = 0;
else
    m_min=2*mc-m_max;
end

m = linspace(m_min/1.1,m_max*1.1,100);
B = olfert.mass2mob(m,mass_mob_pref,mass_mob_exp,P,T);

for i=1:length(m)
    
    [trans(i)]=olfert.CPMA_nodiff_model(m(i),B(i),q,Q,omega,V,CPMA_dims);
    %[trans_D(i)]=olfert.CPMA_diff_model(m(i),B(i),q,Q,omega,V,T,CPMA_dims);
    [trans_Df(i)]=olfert.CPMA_diff_model_fixedgrid(m(i),B(i),q,Q,omega,V,T,CPMA_dims);
    [Omega(i)]=olfert.triangular_transfer_fcn(1/Rm,m(i)/mc);
    
end

figure(3);
%plot(m,trans,m,trans_Df);
% 
% plot(m,trans,m,trans_Df,m,Omega);
hold on;
plot(m,trans_Df);
hold off;

%figure
%plot(m,trans_Df-trans_D);

