function [trans]=CPMA_nodiff_model(m,B,q,Q,omega,V,CPMA_dims)



omega1 = omega/((CPMA_dims.r_hat^2-CPMA_dims.omega_hat)/(CPMA_dims.r_hat^2-1)+CPMA_dims.r1^2*(CPMA_dims.omega_hat-1)/(CPMA_dims.r_hat^2-1)/CPMA_dims.rc^2);


%some controls
data_plot = 0;              % 1 to plot data
delta_r =   1e-9;            %disturbance to particle (1 nm)
num_part=   2;                 %number of particles;


A = pi*(CPMA_dims.r2^2-CPMA_dims.r1^2);       %cross sectional area of APM
v_flow = Q/A;               %average flow velocity



%pick the range of starting points
r_start = linspace(CPMA_dims.r1+delta_r,CPMA_dims.r2-delta_r,num_part);
zspan = [CPMA_dims.L 0];

%for plots...
if data_plot == 1
    figure;
    hold on;
end
for i=1:length(r_start)
    %ODE solver ode23s - Runge-Kutta method (since it is stiff???)
    r_initial = [r_start(i)];
    options = odeset('Events',@olfert.events_CPMS,'RelTol',1e-6);

        [z,r_ode,zE,r_odeE,iE]=ode23s(@olfert.CPMA_ode_fcn,zspan,r_initial,options,v_flow,CPMA_dims.r2,CPMA_dims.r1,omega1,CPMA_dims.omega_hat,B,m,q,V);                          

    if data_plot == 1
        hold on;
        plot(z,r_ode);
        xlabel('Axial Distance (m)');
        ylabel('Radial Direction (m)');
        title('Trajectory of Particle in APM'); 
        ylim([CPMA_dims.r1 CPMA_dims.r2]);
        xlim([0 L]);
    end

    n=length(z);
    z_end(i)=z(n);
    r_end(i)=r_ode(n);
end

del=(CPMA_dims.r2-CPMA_dims.r1)/2;
rho_h=(r_end(2) - CPMA_dims.rc)/del;
rho_l=(r_end(1) - CPMA_dims.rc)/del;
trans=(rho_h-rho_l)*(3-rho_h^2-rho_h*rho_l-rho_l^2)/4;

end

