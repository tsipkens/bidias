

%% CPC counting efficiency
% http://www.sciencedirect.com/science/article/pii/S0021850207000705
% from Wiedensohler ~2000. These parameters are for the 3786 used in the
% 2000 Wiedensohler paper.

a = 0.994;
b = 1.014;
d1 = 4.375E-9; % [nm]
d2 = 0.712E-9; % [nm]

d0 = d2*log(b/(a-1)) + d1;

CPCEff = a - b*((1 + exp((d - d1)./d2)).^(-1));
transform = d > d0;
CPCEff = CPCEff.*transform;


%% Transmission losses
% also assuming mu <0.02. This is the Gormley Kennedy Penetration equation

eta = pi*D*L/Q_a;
Pen = 1.0 - 2.5638*eta.^(2/3) + 1.2*eta + 1.767*eta.^(4/3); % Gormley Kennedy 1949


