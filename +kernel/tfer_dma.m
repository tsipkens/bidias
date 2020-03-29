
% TFER_DMA  Evaluates the transfer function of a differential mobility analyzer.
% Author: Timothy Sipkens, 2018-12-27
% Adapted: Buckley et al. (2017) and Olfert group
% 
% Inputs:
%   d_star          Particle diameter, measurement set point for DMA [m]
%   d               Particle diameter, points in integral, can be vector [m]
%   z               Integer charge state, scalar
%   prop            DMA properties, struct, generated using prop_DMA function (optional)
%   opts.diffusion  Indicates whether to include diffusion, boolean (optional)
%       .solver     Indicates the method by which diffusion is calculated (optional)
%       .param      String indicated which parameter set to use (see prop_DMA.m)
%
% Outputs:
%   Omega           Transfer function
%   Zp_tilde        Non-dimensional electrical mobility, vector
%=========================================================================%

function [Omega,Zp_tilde] = tfer_dma(d_star,d,z,prop,opts)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('opts','var'); opts = []; end % initialize options struct
if ~exist('prop','var'); prop = []; end

if ~isfield(opts,'solver'); opts.solver = 'fullydeveloped'; end
if ~isfield(opts,'diffu'); opts.diffusion = 1; end
if isempty(prop); prop = kernel.prop_dma(opts); end
%-------------------------------------------------------------------------%


%-- Physical constants ---------------------------------------------------%
kb = 1.38064852e-23; % Boltzmann constant [m^2 kg s^-2 K^-1]
e = 1.6022E-19; % electron charge [C]


%-- Evaluate particle mobility -------------------------------------------%
if strcmp(opts.solver,'Buckley')
    [B,Zp] = tfer_pma.dm2zp(d,z); % evaluate electrical mobility (Davies)
    [~,Zp_star] = tfer_pma.dm2zp(d_star);
else
    [B,Zp] = tfer_pma.dm2zp(d,z,prop.T,prop.p); % evaluate electrical mobility (Kim et al.)
    [~,Zp_star] = tfer_pma.dm2zp(d_star,1,prop.T,prop.p);
end
Zp_tilde = Zp./Zp_star; % array of non-dimensional mobilities


%-- Evaluate transfer function -------------------------------------------%
if opts.diffusion
    switch opts.solver
        case 'buckley' % evaluation from Buckley et al.
            D = prop.D(B).*z; % diffusion
            sigma = (prop.G_DMA*2*pi*prop.L.*D./prop.Q_c).^0.5;
            
        case 'plug' % Plug flow, Stolzenburg, 2018
            V = (prop.Q_c/(2*pi*Zp_star*prop.L))*log(prop.R2/prop.R1); % Classifier Voltage (TSI DMA 3080 Manual Equation B-5)
            sigma_star = sqrt(((kb*prop.T)/(z*e*V))*prop.G_DMA); % Stolzenburg Manuscript Equation 20
            sigma = sqrt(sigma_star^2.*Zp_tilde); % Stolzenburg Manuscript Equation 19
            
        case {'olfert','fullydeveloped'} % Olfert laboratory; Fully developed flow, Stolzenburg, 2018
            V = (prop.Q_c/(2*pi.*Zp_star.*prop.L)).*log(prop.R2/prop.R1); % Classifier Voltage (TSI DMA 3080 Manual Equation B-5)
            sigma_star = sqrt(((kb*prop.T)./(z.*e*V))*prop.G_DMA); % Stolzenburg Manuscript Equation 20
            sigma = sqrt(sigma_star^2.*Zp_tilde); % Stolzenburg Manuscript Equation 19
            
        otherwise
            disp('Invalid solver specified.');
            return;
            
    end
    
    epsilon = @(x) x.*erf(x)+exp(-x.^2)./(pi^0.5);
    
    %-- Standard DMA transfer function for diffusion --%
    %	Computes the diffusive transfer function for the DMA
    %	based on Stolzenberg's 1988 analysis. 
    Omega = sigma./(2^0.5*prop.bet*(1-prop.del)).*(...
        epsilon((Zp_tilde-1-prop.bet)./(2^0.5*sigma))+...
        epsilon((Zp_tilde-1+prop.bet)./(2^0.5*sigma))-...
        epsilon((Zp_tilde-1-prop.bet*prop.del)./(2^0.5*sigma))-...
        epsilon((Zp_tilde-1+prop.bet*prop.del)./(2^0.5*sigma)));
    
else % simpler evaluation for the case of exluding diffusion
    Omega = 1/(2*prop.bet*(1-prop.del)).*(abs(Zp_tilde-1+prop.bet)+...
        abs(Zp_tilde-1-prop.bet)-...
        abs(Zp_tilde-1+prop.bet*prop.del)-...
        abs(Zp_tilde-1-prop.bet*prop.del));
end
%-------------------------------------------------------------------------%


%-- Remove negative and small values --%
negTransform = Omega > 0;
Omega = negTransform.*Omega;
smallTransform = Omega > 1e-14;
Omega = smallTransform.*Omega;


Omega = Omega'; % transpose data to output desired format

end

