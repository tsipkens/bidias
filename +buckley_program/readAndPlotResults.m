function [ChiSquared,dn_ddpdmp,dM_ddp,dn_ddp,dp,mp] = readAndPlotResults(dataName)
cd('Results')

ChiSquared = csvread([dataName,'ChiSquared.csv']);
dn_ddpdmp = csvread([dataName,'distribution.csv']);
dM_ddp = csvread([dataName,'dM_ddp.csv']);
dn_ddp = csvread([dataName,'dn_ddp.csv']);
dp = csvread([dataName,'dp.csv']);
mp = csvread([dataName,'mp.csv']);

cd ../

dp = dp.*1E9; % convert diameters to nm
mpDa = mp.*6.022E26; % convert mass from kg to Da

maxZ = max(max(dn_ddpdmp));
Zlarge = dn_ddpdmp > 0.000001*maxZ;
Z = dn_ddpdmp.*Zlarge;


%% Plot dn/ddpdmp
contourf(dp,mpDa,log10(transpose(Z)),30,'LineStyle','none')
xlabel('d_p [nm]')
ylabel('m_p [Da]')
set(gca,'fontsize',18)
% axis([0 450 0 0.75E-17]) % can set axis limits if desired with this
%                          % command

%% Plot dM/ddp
dM_ddpSmoothed = simpleSmooth(dM_ddp,100);
figure
plot(dp,dM_ddpSmoothed)
xlabel('d_p [nm]')
ylabel('dM/dd_p [kg/cm^{-3}]')


%% Plot dn/ddp
dn_ddpSmoothed = simpleSmooth(dn_ddp,100);
figure
plot(dp, dn_ddpSmoothed)
xlabel('d_p [nm]')
ylabel('dn/dd_p [cm^{-3}]')

%% Plot Chi-Squared Error

% figure
% plot(ChiSquared)
% xlabel('k (Twomey iterations)')
% ylabel('Chi-Squared')
