function [B,d_m] = mass2mob(m,mass_mob_pref,mass_mob_exp,P,T)

for i=1:length(m)
    d_m(i) = (m(i)/mass_mob_pref)^(1/mass_mob_exp);
    B(i) = olfert.dp2mobility_v3(d_m(i),P,T);
end