
%-- In-house, demonstration phantom --------------%
phantom(1).rho = 1500; % 19300; % density of gold-ish
phantom(1).Dm = Dm1_vec(ii);
phantom(1).dg = 120;
phantom(1).sg = 1.6;
phantom(1).k = (phantom(1).rho*pi/6).*...
    phantom(1).dg^(3-phantom(1).Dm);
phantom(1).sm = 1.4;
phantom(1).opt_m = 'logn';

% phantom(2).rho = 500;
% phantom(2).Dm = 2.3;
% phantom(2).dg = 200;
% phantom(2).sg = 1.4;
% phantom(2).k = (phantom(2).rho*pi/6).*...
%     phantom(2).dg^(3-phantom(2).Dm);
% phantom(2).sm = 1.3;
% phantom(2).opt_m = 'logn';



