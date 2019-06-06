function [CPMA_dims] = CPMA_Camb_dims(omega_hat)

CPMA_dims.gap     =   .001;               %distance between cylinders (m)
CPMA_dims.r1      =   .060;               %inner radius of APM (m)
CPMA_dims.r2      =   CPMA_dims.r1+CPMA_dims.gap;             %outer radius of APM (m)
CPMA_dims.r_hat   =   CPMA_dims.r1/CPMA_dims.r2;
CPMA_dims.omega_hat   =   omega_hat;
CPMA_dims.rc      =   (CPMA_dims.r2+CPMA_dims.r1)/2;          %the center radius
CPMA_dims.L       =   .200;               %length of electrode (cylinder) (m)
CPMA_dims.k       =   1.3806503*10^-23;   %boltzmann's constant