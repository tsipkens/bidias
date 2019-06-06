function [T_old]=CPMA_diff_model(m,B,q,Q,omega,V,T,CPMA_dims)

omega1 = omega/((CPMA_dims.r_hat^2-CPMA_dims.omega_hat)/(CPMA_dims.r_hat^2-1)+CPMA_dims.r1^2*(CPMA_dims.omega_hat-1)/(CPMA_dims.r_hat^2-1)/CPMA_dims.rc^2);



r_size  =   100;                %number of r steps
plot_data=  0;

%%% FIND PARTICLE DIMENSIONS %%%

D = 1.3806488e-23*T*B;                                  %Diffusion coefficient



%%%% find the forces %%%%
r=linspace(CPMA_dims.r1,CPMA_dims.r2,r_size);
%calculate the velocity profile for the couette flow
NU = CPMA_dims.r1/CPMA_dims.rc;
MU = omega/omega1;
v = omega1 * (NU^2 - MU) / (NU^2 - 1) .* r + omega1 * CPMA_dims.r1^2 * (MU - 1) / (NU^2 - 1) ./ r;


Fc = m .* v.^2 ./ r;       %centrifugal force (N)
Fe = q .* V ./ r ./ log(CPMA_dims.r2 / CPMA_dims.r1);    %the electrostatic force (N)
delta_F = Fc-Fe;

%transform the cylinder geometery to flat plates dimensions
y = r-CPMA_dims.rc;

X = [ones(size(y'))  y'  y'.^2];
a = X\delta_F';

if plot_data == 1
    %figure;
    %plot(y,Fc,y,Fe);
    %legend('Centripetal Force','Electrostatic Force')
    
    %figure;
    %plot(y,v)
    
    figure
    plot(r,delta_F,'.')
    %hold on;
    %ff = (-h:0.0001:h)';
    %Y = [ones(size(ff))  ff  ff.^2]*a;
    %plot(ff,Y)
end


%need flow rate
A = pi*(CPMA_dims.r2^2-CPMA_dims.r1^2);     %cross sectional area of APM
v_flow = Q/A;           %average flow velocity


%calculate the non-dimensional variables
a2_tld = (CPMA_dims.gap/2)^3 *B / D * a(3);
a1_tld = (CPMA_dims.gap/2)^2 *B / D * a(2);
a0_tld = (CPMA_dims.gap/2)   *B / D * a(1);
x_max = 2/3 * D * CPMA_dims.L / v_flow / (CPMA_dims.gap/2)^2;


%the initial number of y & x steps
m=50;
H=2*m+1;
K=H;

%some constants
tol = .1/100;    %tolerance
error = 1;      %initialize error
j=1;            %initialize run counter
T_old= 1;
loop_control = 0;


while error>tol
    
    h=2/(H-1);      %step size in y~
    k=x_max/(K-1);  %step size in x~
    
    %starting vector (initial conditions)
    n_old = linspace(1,1,H)';
    
    
    %calculate a velocity vector for the parabolic flow
    vel_flow(H)=0;      %preallocate for speed
    for i = 1:H
        y = -1 + (i-1)*h;
        vel_flow(i) =  (1 - y^2); %I don't need the 3/2, it gets divided out
    end
    
    n=zeros(H,1);
    count = 1;
    n(:,count)=n_old;
    
    
    %calculate the initial concentration of particles
    
    temp = n_old.*vel_flow';
    flux_in = trapz(temp);
    
    
    
    x=0;
    while x <= x_max
        x=x+k;
        
        %%%%set up the systems of equations to solve%%%%
        for i = 2:H-1
            y = -1 + (i-1)*h;
            W = 2*a2_tld*y + a1_tld;
            Q = a2_tld*y^2 + a1_tld*y + a0_tld;
            U = 1 - y^2;
            
            b(i) = U/k + 1/h^2 + W/2;
            a(i-1) = -1/2/h^2 - Q/4/h;
            c(i) = Q/4/h - 1/2/h^2;
            RHS(i) = n_old(i-1)*(1/2/h^2 + Q/4/h) + n_old(i)*(U/k - 1/h^2 - W/2) + n_old(i+1)*(1/2/h^2 - Q/4/h);
            
        end
        
        %%%%% the boundary conditions %%%%%
        %for concentration of zero at the wall
        b(1) = 1;
        b(H) = 1;
        c(1) = 0;
        a(H-1) = 0;
        RHS(1) = 0;
        RHS(H) = 0;
        
        %this is a tridiagonal matrix solver
        n_new = olfert.tridiag1(a,b,c,RHS);
        
        %Calculate the transfer function
        temp2 = n_new.*vel_flow;
        flux_out = trapz(temp2);
        T_new=flux_out/flux_in;
        
        %loop control
        count = count+1;
        n(:,count)=n_new';
        
        
        if T_new < 0.0005
            loop_control =1;
            break
        end
        
        n_old = n_new;
    end
    
  
    if loop_control ==1;
      T_old=0;
       break
    end
    
    error = max(abs(T_old-T_new));
    T_old=T_new;
    j=j+1;
    m = 2*m;
    
    H=2*m+1;
    K=H;

end

if plot_data == 1
    [R,C]=size(n);
    plt_msh=50;
    xx=linspace(0,x,C);
    yy=linspace(-1,1,R);
    xi=linspace(0,x,plt_msh);
    yi=linspace(-1,1,plt_msh);
    X=zeros(R,C);
    Y=zeros(R,C);
    XI=zeros(plt_msh,plt_msh);
    YI=zeros(plt_msh,plt_msh);
    for i=1:C
        X(:,i)=xx(i);
    end
    for i=1:R
        Y(i,:)=yy(i);
    end
    for i=1:plt_msh
        XI(:,i)=xi(i);
    end
    for i=1:plt_msh
        YI(i,:)=yi(i);
    end
    nI = interp2(X,Y,n,XI,YI);
end

if plot_data == 1
    figure;
    mesh(XI,YI,nI)
    xlabel('x^~');
    ylabel('y^~');
    zlabel('n^~');
    xlim([0 x_max]);
    view([45 45]);
    title('Particle Concentration in APM');
    %zlim([0 1])
end
