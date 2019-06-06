%This m-file plots the transfer function given by Stolzenburg and McMurry
%(2008) AST vol 42, page 421-432.

%function DMA_D_Stolzenburg()

Q_s = 0.3/60/1000; %Sample flow, m^3/s
Q_a = 0.3/60/1000; %Aerosol flow, m^3/s
Q_sh = 3/60/1000; %Sheath flow, m^3/s
Q_exh = 3/60/1000; %Exhaust flow, m^3/s
L=0.44369; %length, m
R1=0.00937; %inner radius, m
R2=0.01961; %outer radius, m
T=298; %temperature, K
k = 1.3806503*10^-23; %Boltzmann's constant
q = 1.60217646*10^-19; %elementary charge
dm_star = d_star;%100e-9; % Mobility Diameter Setpoint in m %%%%%%%%%

beta = (Q_s + Q_a)/(Q_sh + Q_exh);
delta = (Q_s - Q_a)/(Q_s + Q_a);

Zp_tilde = linspace(.6, 1.4, 101);%kernel.dp2elecmob(d.*1e-9,1,1,T);%linspace(.6, 1.4, 101);
Zp_star = olfert.dp2elecmob(dm_star,1,1,T);
Zp = Zp_tilde*Zp_star;
B = Zp/q;
D = k*T*B;

Omega_ND = 1/(2*beta*(1-delta))*(abs(Zp_tilde-(1+beta)) + ...
    abs(Zp_tilde-(1-beta)) - abs(Zp_tilde-(1+beta*delta)) - abs(Zp_tilde-(1-beta*delta)));

figure(1);
plot(Zp_tilde,Omega_ND);
xlabel('Z_p/Z_p*');
ylabel('\Omega');
ylim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate the diffusion factor %%%%%%%%%%%%
gamma = (R1/R2)^2; 
%%%%%kappa = L*R2/(R2^2-R1^2);
kappa = L / R2; % Stolzenburg Manuscript Eqn 9

f=@(omega)((1-omega)^2*log(gamma)/2/(1-gamma)+omega*log(omega)+(1-omega))/((1+gamma)*log(gamma)/2+(1-gamma))-beta*(1-delta)/2/(1+beta);
guess=0.9;
omega_a = fzero(f,guess);

f=@(omega)(1-((1-omega)^2*log(gamma)/2/(1-gamma)+omega*log(omega)+(1-omega))/((1+gamma)*log(gamma)/2+(1-gamma)))-beta*(1+delta)/2/(1+beta);
guess=0.3;
omega_s = fzero(f,guess);

A = (-.5*(1+gamma)*log(gamma)-(1-gamma))^(-1);
I_omega_a = A^2/(1-gamma)*(-omega_a^2*((1-gamma)*log(omega_a)-(1-omega_a)*log(gamma))^2/2+(omega_a^2*(1-gamma)/2+omega_a^3*log(gamma)/3)*((1-gamma)*log(omega_a)-(1-omega_a)*log(gamma))+(1-omega_a^2)*(1-gamma)^2/4+5*(1-omega_a^3)*(1-gamma)*log(gamma)/18+(1-omega_a^4)*(log(gamma))^2/12);
I_omega_s = A^2/(1-gamma)*(-omega_s^2*((1-gamma)*log(omega_s)-(1-omega_s)*log(gamma))^2/2+(omega_s^2*(1-gamma)/2+omega_s^3*log(gamma)/3)*((1-gamma)*log(omega_s)-(1-omega_s)*log(gamma))+(1-omega_s^2)*(1-gamma)^2/4+5*(1-omega_s^3)*(1-gamma)*log(gamma)/18+(1-omega_s^4)*(log(gamma))^2/12);

%%%%%G_DMA = (4*(1+beta)^2*(I_omega_s-I_omega_a)+(omega_a-omega_s)/kappa^2)/(1-gamma); % this is not exactly the same as what is found in Stolzenburg thesis but it gives a similar answer
%%%%%D_tilde = 4*pi*L*D/(Q_sh+Q_exh);
%%%%%sigma = sqrt(G_DMA*D_tilde);
G_o_corr = (4 * (1 + beta)^2) * (I_omega_s - I_omega_a) / (1 - gamma) + (omega_a - omega_s) / kappa^2; % Stolzenburg Manuscript Equation 8 and 12
G_DMA = G_o_corr * log(R2 / R1); % Stolzenburg Manuscript Equation 23
V = (Q_sh / (2 * pi * Zp_star * L)) * log(R2 / R1); % Classifier Voltage (TSI DMA 3080 Manual Equation B-5)
sigma_star = sqrt(((k * T) / (1 * q * V)) * G_DMA); % Stolzenburg Manuscript Equation 20
sigma = sqrt(sigma_star^2 .* Zp_tilde); % Stolzenburg Manuscript Equation 19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Zp_tilde)
Omega_D(i) = sigma(i)/(sqrt(2)*beta*(1-delta))*(epsilon((Zp_tilde(i)-(1+beta))/sqrt(2)/sigma(i)) + ...
    epsilon((Zp_tilde(i)-(1-beta))/sqrt(2)/sigma(i)) - epsilon((Zp_tilde(i)-(1+beta*delta))/sqrt(2)/sigma(i)) -...
    epsilon((Zp_tilde(i)-(1-beta*delta))/sqrt(2)/sigma(i)));
end

hold on;
plot(Zp_tilde,Omega_D);
hold off;

%end

function [epln]=epsilon(x)
epln = x*erf(x)+exp(-x^2)/sqrt(pi);
end
