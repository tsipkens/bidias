function [G] = DMATransfer(invZpArray, invZpStar, dp, n)
% function [DMATfer] = KernelCalc_InvMobility(invZpArray, Qsh, ZpInput)

%% DMA Parameters (USER INPUT)

Qsh = 4.89E-3/60; % m^3/sec
Qae = 1.02E-3/60; % m^3/sec
L = 0.4444; % [m]
R2 = 0.01961; % [m]
R1 = 0.00937; % [m]
Ltube = 0.1; % length of tube, change only if accounting for diffusional losses [m]

%%%%% Initial Parameters:

kb = 1.38064852E-23; % boltzmann [m^2 kg s^-2 K^-1]
T = 294; % temp [K]
e = 1.6022E-19; % electron charge [C]

%%%%% Mobility calculation:

invZpArray = invZpArray./n; % account for charge state
ZpStar = 1/invZpStar;

%%%%% Transfer Function Calculation:

y = (R1/R2)^2;
kap = L*R2/(R2^2 - R1^2);
Iy = (0.25*(1-y^2)*(1-y)^2 + (5/18)*(1 - y^3)*(1 - y)*log(y) + (1/12)*(1-y^4)*(log(y))^2)/((1-y)*(-0.5*(1+y)*log(y) - (1-y))^2);
B = (Qae/Qsh);
Z_p = (invZpArray.*ZpStar).^(-1); % array of non-dimensional mobilities
GDMA = 4*(1 + B)^2*(Iy + (2*(1+B)*kap)^(-2))/(1-y);
D = kb*T./(invZpArray*e);
sig = (GDMA*2*pi*L*D/Qsh).^0.5;

DMATfer = (2^0.5*B)^(-1)*sig.*(Old.eps1((Z_p - (1 + B))./(2^0.5*sig)) + Old.eps1((Z_p - (1 - B))./(2^0.5*sig)) - 2.*Old.eps1((Z_p - 1)./(2^0.5*sig)));
% The line above computes the diffusive transfer function for the DMA
% based on Stolzenberg's 1988 analysis for previously defined Zp range


%% CPC Counting Efficiency
% http://www.sciencedirect.com/science/article/pii/S0021850207000705
% from Wiedensohler ~2000. These parameters are for the 3786 used in the 2000
% Wiedensohler paper.

a = 0.994;
b = 1.014;
d1 = 4.375E-9; % [nm]
d2 = 0.712E-9; % [nm]

d0 = d2*log(b/(a-1)) + d1;

CPCEff = a - b*((1 + exp((dp - d1)./d2)).^(-1));
transform = dp > d0;
CPCEff = CPCEff.*transform;

%% Transmission Losses
% also assuming mu <0.02. This is the Gormley Kennedy Penetration equation

eta = pi*D*Ltube/Qae;
Pen = 1.0 - 2.5638*eta.^(2/3) + 1.2*eta + 1.767*eta.^(4/3); % Gormley Kennedy 1949

%% Kernel Calculation

% G = DMATfer.*Pen.*CPCEff; % Uncomment this function to include
% penetration and cpc counting efficiency
G = DMATfer; % comment this function to include penetration and cpc counting efficiency
G = transpose(G);

negTransform = G > 0;
G = negTransform.*G;
smallTransform = G > 1E-14;
G = smallTransform.*G;