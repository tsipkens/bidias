
% TFER_DMA  Evaluates the transfer function of a differential mobility analyzer.
% 
%  OMEGA = kernel.tfer_dma(D_STAR,D,Z) uses the mobility diameter set points
%  specified by D_STAR [m] and evalautes the DMA transfer function at D
%  [m] for an integer charge state of Z (a scalar, integer). Uses default
%  properties specified by kernel.prop_dma. Explicitly stating prop_dma is
%  preferred. Output is transfer function, OMEGA. 
%  
%  OMEGA = kernel.tfer_dma(D_STAR,D,Z,PROP_DMA) explicitly specified the
%  properties of the DMA (e.g., the radii) as a data structure. This is the
%  preferred usage to the previous call. For structure of PROP_DMA, see
%  kernel.prop_dma(...). 
%  
%  OMEGA = kernel.tfer_dma(...,OPTS) adds an options structure with the
%  field specified below: 
%   OPTS.diffusion  Indicates whether to include diffusion, boolean (optional)
%       .solver     Indicates the method by which diffusion is calculated (optional)
%       .param      String indicated which parameter set to use (see prop_DMA.m)
%  
%  [OMEGA,ZP_TILDE] = kernel.tfer_dma(...) adds an output contianing the
%  non-dimensional electrical mobility as a vector.
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens, 2018-12-27
%  ADAPTED FROM: Buckley et al. (2017) and Olfert group


function [Omega,Zp_tilde] = tfer_dma(d_star,d,z,prop,opts)


%-- Parse inputs ---------------------------------------------------------%
addpath tfer_pma; % add mat-tfer-pma package to MATLAB path

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
    [B, Zp] = dm2zp(d,z); % evaluate electrical mobility (Davies)
    [~, Zp_star] = dm2zp(d_star);
else
    [B, Zp] = dm2zp(d, z, prop.T, prop.p); % evaluate electrical mobility (Kim et al.)
    [~, Zp_star] = dm2zp(d_star, 1, prop.T, prop.p);
end
Zp_tilde = Zp ./ Zp_star; % array of non-dimensional mobilities


%-- Evaluate transfer function -------------------------------------------%
if opts.diffusion
    switch opts.solver
        case 'buckley' % evaluation from Buckley et al.
            D = prop.D(B) .* z; % diffusion
            sigma = (prop.G_DMA*2*pi*prop.L .* D ./ prop.Q_c) .^ 0.5;
            
        case 'plug' % Plug flow, Stolzenburg, 2018
            V = (prop.Q_c ./ (2*pi .* Zp_star .* prop.L)) .* log(prop.R2 / prop.R1); % Classifier Voltage (TSI DMA 3080 Manual Equation B-5)
            sigma_star = sqrt(((kb*prop.T) / (z*e*V)) * prop.G_DMA); % Stolzenburg Manuscript Equation 20
            sigma = sqrt(sigma_star^2.*Zp_tilde); % Stolzenburg Manuscript Equation 19
            
        case {'olfert', 'fullydeveloped'} % Olfert laboratory; Fully developed flow, Stolzenburg, 2018
            V = (prop.Q_c ./ (2*pi .* Zp_star .* prop.L)) .* log(prop.R2 / prop.R1); % Classifier Voltage (TSI DMA 3080 Manual Equation B-5)
            sigma_star = sqrt(((kb*prop.T) ./ (z .* e .* V)) .* prop.G_DMA); % Stolzenburg Manuscript Equation 20
            sigma = sqrt(sigma_star .^ 2 .* Zp_tilde); % Stolzenburg Manuscript Equation 19
            
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

