
% GEN_1D_DMA  Evaluates the transfer function of a differential mobility analyzer.
% Author:	Timothy Sipkens, 2020-03-09
% 
% Inputs:
%   d_star              Particle diameter, measurement set point for DMA [m]
%   d                   Particle diameter, points in integral, can be vector [m]
% 
%   varargin{1} = 
%       prop	DMA properties, struct, generated using prop_DMA function (optional)
%   
%   varargin{2} = 
%       opts.diffusion  Indicates whether to include diffusion, boolean (optional)
%           .solver     Indicates the method by which diffusion is calculated (optional)
%           .param      String indicated which parameter set to use (see prop_DMA.m)
%
% Outputs:
%   Omega           Transfer function
%   Zp_tilde        Non-dimensional electrical mobility, vector
%=========================================================================%

function [Omega] = gen_1d_dma(d_star,d,varargin)

n_b = length(d_star);
n_i = length(d);

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';
f_z = sparse(kernel.tfer_charge(d.*1e-9,z_vec)); % get fraction charged for d
n_z = length(z_vec);


%== Evaluate DMA transfer function =======================================%
tools.textheader('Computing DMA kernel');
Omega = sparse(n_b,n_i);
for kk=1:n_z
    t0 = zeros(n_b,n_i); % pre-allocate for speed
    
    for ii=1:n_b % loop over d_star
        t0(ii,:) = kernel.tfer_dma(...
            d_star(ii).*1e-9,...
            d.*1e-9,...
            z_vec(kk),varargin{:});
    end
    
    t0(t0<(1e-7.*max(max(t0)))) = 0;
        % remove numerical noise in kernel
    
    Omega = Omega+f_z(kk,:).*sparse(t0);
end
tools.textheader();


end
