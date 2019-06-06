
%-- In-house, demonstration phantom --------------%
%-{
phantom(1).rho = 12000; % 19300; % density of gold-ish
phantom(1).Dm = 3;
phantom(1).dg = 50;
phantom(1).sg = 1.4;
phantom(1).k = (phantom(1).rho*pi/6).*...
    phantom(1).dg^(3-phantom(1).Dm);
phantom(1).sm = 1.3;
phantom(1).opt_m = 'logn';

phantom(2).rho = 500;
phantom(2).Dm = 2.3;
phantom(2).dg = 200;
phantom(2).sg = 1.4;
phantom(2).k = (phantom(2).rho*pi/6).*...
    phantom(2).dg^(3-phantom(2).Dm);
phantom(2).sm = 1.3;
phantom(2).opt_m = 'logn';
%}


%-- Buckley/Hogan phantom ----------------------%
%{
phantom(1).rho = 10000;
phantom(1).Dm = 3;
phantom(1).dg = 200;
phantom(1).sg = 1.5;
phantom(1).k = (phantom(1).rho*pi/6).*...
    phantom(1).dg^(3-phantom(1).Dm);
phantom(1).sm = 0.15;
phantom(1).opt_m = 'norm';

phantom(2).rho = 1000;
phantom(2).Dm = 3;
phantom(2).dg = 300;
phantom(2).sg = 2.2;
phantom(2).k = (phantom(2).rho*pi/6).*...
    phantom(2).dg^(3-phantom(2).Dm);
phantom(2).sm = 0.15;
phantom(2).opt_m = 'norm';
% Buckley plot --> 1e9 Da = 1.66054 fg
               % --> 0.415 fg to 6.642 fg
               % --> 10^-0.3820 fg to 10^0.8223 fg
%}


%-- Phantom based on Olfert data ---------%
%{
phantom(1).Dm = 2.3;
phantom(1).dg = 125;
phantom(1).sg = 1.6;
phantom(1).k = 9400;
phantom(1).rho = 6*phantom(1).k/pi*phantom(1).dg^(phantom(1).Dm-3);
phantom(1).sm = 1.5;
phantom(1).opt_m = 'logn';
%}


%-- Narrow phantom -----------------%
%{
phantom.rho = 1000;
phantom.Dm = 3; % 2.3
phantom.dg = 125;
phantom.sg = 1.5;
phantom.k = (phantom.rho*pi/6).*...
    phantom.dg^(3-phantom.Dm);
% phantom.sm = 1.02;
phantom.sm = 1.05;
phantom.mg = 1e-9*phantom.k*phantom.dg^phantom.Dm;
phantom.opt_m = 'logn';
%}


