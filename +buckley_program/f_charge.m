function [f0, f1, f2] = f_charge(dp)
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


%% f0 nonconducting particles
% a0 = -1.212;
% a1 = 1.1068;
% a2 = -0.2934;
% a3 = 0.0169;

%% f0 conducting particles w/ image potential
a0 = -0.3880;
a1 = 0.4545;
a2 = -0.1634;
a3 = 0.0091;


a_i = [a0,a1,a2,a3];

exponent = 0;
for i = 1:4
    exponent = exponent + a_i(i)*(log(dp).^(i-1));
end

f0 = exp(exponent);
              

%% f1

% nonconducting
% a0 = -16.704;
% a1 = 7.5438;
% a2 = -1.1938;
% a3 = 0.0589;

% conducting
a0 = -8.0157;
a1 = 3.2536;
a2 = -0.5018;
a3 = 0.0223;

a_i = [a0,a1,a2,a3];

exponent = 0;
for i = 1:4
    exponent = exponent + a_i(i)*(log(dp).^(i-1));
end

f1 = exp(exponent);

%% f2

clear [a0, a1,a2,a3,a_i, exponent]

% a0 = -71.051;
% a1 = 31.209;
% a2 = -4.6696;
% a3 = 0.2301;

% conducting particles

a0 = -40.714;
a1 = 17.487;
a2 = -2.6146;
a3 = 0.1282;

a_i = [a0,a1,a2,a3];

exponent = 0;
for i = 1:4
    exponent = exponent + a_i(i)*(log(dp).^(i-1));
end

f2 = exp(exponent);

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