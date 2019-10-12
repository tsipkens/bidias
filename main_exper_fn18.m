
clear;
clc;
close all;


%=========================================================================%
%-- Load colour schemes --------------------------------------------------%
addpath('cmap');
cm = load_cmap('YlGnBu',255);
cm_alt = cm;
load('inferno.mat');
cm = cm(40:end,:);
% load('matter.mat');
cm_b = cm;
load('viridis.mat');


files = {'20180528_A_L9.mat',...
    '20180528_G_M9.mat',...
    '20180601_E_M9.mat',...
    '20180528_F_NS_M9.mat',...
    '20180529_C_AB_M9.mat',...
    '20180528_B_H9.mat',...
    '20180528_E_RU-KM2.mat',...
    '20180528_D_RU-KM1.mat',...
    '20180529_D_EC_AC.mat',...
    '20180529_B_BK_WO.mat',...
    '20180528_C_BK-BR.mat',...
    '20180529_E2_EC_AC27.mat'};
fuel = {'L9','M9','M9','NS_M9','AB_M9',...
    'H9','RU_KM2','RU_KM1','EC_AC','BK_WO',...
    'BK_BR','EC_AC27'};
% NOTE: EC-AS1 (20180601_A_EC_AS1.mat) has anomalous data
data0(length(files)) = struct;
data2(length(files)) = struct;
data3(length(files)) = struct;

for ff=1:length(files)
    
    disp(' ');
    disp(['<== PROCESSING FUEL ',num2str(ff),' OF ',...
        num2str(length(files)),' (',fuel{ff},') =====================>']);
    
    %=========================================================================%
    %-- Load experimental data -----------------------------------------------%
    load(['..\data\Soot Data FlareNet 18\',files{ff}]);
    
    %-- Reformat data --------------------------------------------------------%
    data = data';
    b_max = max(max(data));
    b = data(:)./b_max;
    edges_b = {data_m,data_d'};
    grid_b = Grid(edges_b,...
        [],'logarithmic');
    sig = sqrt(data(:))+max(sqrt(data(:))).*0.01; % estimate of noise
    Lb = diag(sqrt(1./sig));

    figure(5);
    colormap(gcf,cm_b);
    [~,b_m] = grid_b.plot2d_marg(b);

    figure(20);
    n2 = floor(grid_b.ne(1));
    n3 = floor(length(cm_b(:,1))/n2);
    cm_b_mod = cm_b(10:n3:end,:);
    set(gca,'ColorOrder',cm_b_mod,'NextPlot','replacechildren');
    b_plot_rs = reshape(b,grid_b.ne);
    semilogx(grid_b.edges{2},b_plot_rs.*b_max);
    % hold on;
    % semilogx(grid_b.edges{2},sum(b_plot_rs).*Ntot,'k');
    % hold off;


    %=========================================================================%
    %-- Generate A and grid_x ------------------------------------------------%
    ne_x = [50,64]; % number of elements per dimension in x
        % [20,32]; % used for plotting projections of basis functions
        % [40,64]; % used in evaluating previous versions of regularization

    span = [10^-2,50;10,10^3];
    grid_x = Grid(span,...
        ne_x,'logarithmic');

    r_x = grid_x.nodes;
    edges_x = grid_x.edges;
    n_x = grid_x.ne;
    
    disp('Evaluate kernel...');
    prop_PMA = kernel.prop_PMA('FlareNet18');
    A = kernel.gen_A(grid_b,grid_x,prop_PMA,'Rm',5); % generate A matrix based on grid for x and b


    %=========================================================================%
    %-- Perform exponential, rotated regularization --------------------------%
    %-{
    s1 = 1.0;
    s2 = 0.1;
    dtot = @(d1,d2) sqrt(exp(d1).^2+exp(d2).^2);
    theta = -atan2(1,2.5);
    Lex = diag([1/s1,1/s2])*...
        [cos(theta),-sin(theta);sin(theta),cos(theta)];
    lambda_expRot = 1e-3; % 5e-4

    disp('Performing rotated exponential distance regularization...');
    [x_expRot,L] = invert.exp_dist(...
        Lb*A,Lb*b,grid_x.elements(:,2),grid_x.elements(:,1),...
        lambda_expRot,Lex);
    disp('Inversion complete.');
    disp(' ');
    %}
    
    %-- Save copies of data --------------------------------------------------%
    [dr,dr1,dr2] = grid_x.dr;
    dr = dr(1,1); dr1 = dr1(1,1); dr2 = dr2(1,1);
    x_expRot = x_expRot.*b_max;
    x_plot = x_expRot;
    dV_fact = 0.3*16.666667*90;
    data0(ff).fuel = fuel{ff};
    data0(ff).x = x_expRot;
    data0(ff).b = b;
    data0(ff).grid_x = grid_x;
    data0(ff).grid_b = grid_b;
    
    
    %=========================================================================%
    %-- Estimate mass-mobility relation --------------------------------------%
    figure(40);
    colormap(gcf,cm_alt);
    grid_x.plot2d_marg(x_plot.*dr.*dV_fact);
    [data0(ff).Dm,data0(ff).k,data0(ff).rho_100] = ...
        grid_x.fit_mass_mob(x_plot); % line maximzing integration of distr.
    tools.fit_mm_simple;
    tools.fit_mm_dist;
    hold on;
    plot(log10(grid_x.edges{2}),...
        log10(1e-9*100^3*pi/6*510.*(grid_x.edges{2}./100).^2.48),'y');
    hold off;
    xlabel('log_{10}(d)');
    ylabel('log_{10}(m)');
    
    data1(ff) = pha.p;
    data2(ff).m_100 = m_100;
    data2(ff).rho_100 = rho_100;
    data2(ff).Dm = Dm;
    
    %-- Plots for effective density ------------------------------------------%
    [y,grid_rho] = ...
        tools.mass2rho(x_plot.*dr.*dV_fact,grid_x);

    figure(31);
    colormap(gcf,cm_alt);
    grid_rho.plot2d_marg(y);
    xlabel('log_{10}(d)');
    ylabel('log_{10}(\rho_{eff})');
    
    tools.fit_rhod_simple;
    % [m,b] = grid_rho.fit_mass_mob(y,[2,2.8],-0.6);
    grid_rho.plot_line_overlay(...
        [0,log10(6*data0(ff).k/pi)+9],data0(ff).Dm-3,'w');
    rho = 2000; % density of base material
    dpe = 10.^((log10(6*data0(ff).k/(pi*rho))+9)/(3-data0(ff).Dm));
    
    
    figure(32);
    semilogx(d_max,rho_max,'o');
    hold on;
    semilogx(grid_x.edges{2},...
        rho_100.*(grid_x.edges{2}./100).^(Dm-3),'k');
    semilogx(grid_x.edges{2},...
        510.*(grid_x.edges{2}./100).^-0.52,'y');
    hold off;
    
    %-- Plot fit phantom ---%
    figure(41);
    colormap(gcf,cm_alt);
    [~,x_plot_m] = grid_x.plot2d_marg(pha.x.*dr.*dV_fact);
    tools.fit_mm_dist;
    xlabel('log_{10}(d)');
    ylabel('log_{10}(m)');
    
    
    %=====================================================================%
    %-- Traditional mass-mobility fitting ---%
    p_max = polyfit(log10(d_max),log10(grid_b.edges{1}),1);
    figure(60);
    grid_b.plot2d(b_plot_rs);
    colormap(cm_b);
    hold on;
    plot(log10(d_max),log10(grid_b.edges{1}),'wo');
    plot(log10(d_max),polyval(p_max,log10(d_max)),'w');
    hold off;
    
    data3(ff).rho_100 = rho_100;
    data3(ff).Dm = Dm;
    
end

[data1.fuel] = fuel{:};
[data2.fuel] = fuel{:};
[data3.fuel] = fuel{:};



