function [f0, f1, f2] = f_charge_W(dp)
% given an input array of sizes (dp) in m (IMPORTANT), f_charge calculates the fraction of 
% charged particles after having been passed through a bipolar diffusion
% charger. f1 and f2 are the respective outputted charge fractions for z =
% +1 and +2. This calculation is done using the equations and values found
% in Gopalakrishnan et al, 2013. Also, we assume the particles are both
% spherical and non-conducting. hyperlink for Gopalakrishnan paper is
% below:
% http://ac.els-cdn.com/S0021850213001079/1-s2.0-S0021850213001079-main.pdf?_tid=77ac684c-efae-11e6-a2bf-00000aab0f6b&acdnat=1486744509_d43fef7d2ec9a8118569636b58e72e5e

dp = dp.*1E9; % IMPORTANT UNITS: the equations below function if dp is in nm, but this
              % function is designed to take in data in units of meters
              % (since they are often converted from 1/K or 1/Zp, etc).


%% f0
a0 = -0.0003;
a1 = -0.1014;
a2 = 0.3073;
a3 = -0.3372;
a4 = 0.1023;
a5 = -0.0105;

a_i = [a0,a1,a2,a3,a4,a5];

exponent = 0;
for i = 1:6
    exponent = exponent + a_i(i)*(log10(dp).^(i-1));
end

f0 = 10.^(exponent);
              

%% f1

a0 = -2.3484;
a1 = 0.6044;
a2 = 0.4800;
a3 = 0.0013;
a4 = -0.1544;
a5 = 0.0320;

a_i = [a0,a1,a2,a3,a4,a5];

exponent = 0;
for i = 1:6
    exponent = exponent + a_i(i)*(log10(dp).^(i-1));
end

f1 = 10.^(exponent);

%% f2

clear [a0, a1,a2,a3,a_i, exponent]

a0 = -44.4756;
a1 = 79.3772;
a2 = -62.8900;
a3 = 26.4492;
a4 = -5.7480;
a5 = 0.5059;

a_i = [a0,a1,a2,a3,a4,a5];

exponent = 0;
for i = 1:6
    exponent = exponent + a_i(i)*(log10(dp).^(i-1));
end

f2 = 10.^(exponent);

%% f2 (conducting)

% clear [a0, a1,a2,a3,a_i, exponent]
% 
% a0 = -40.714;
% a1 = 17.487;
% a2 = -2.6146;
% a3 = 0.1282;
% 
% a_i = [a0,a1,a2,a3];
% 
% exponent = 0;
% for i = 1:4
%     exponent = exponent + a_i(i)*(log(dp).^(i-1));
% end
% 
% f2 = exp(exponent);