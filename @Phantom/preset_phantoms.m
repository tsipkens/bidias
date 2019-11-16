
% PRESET_PHANTOMS Parameter sets for preset/sample phantoms.
% Author: Timothy Sipkens, 2019-10-31
%=========================================================================%

function [p,type,name] = preset_phantoms(name)

switch name
    case {'demonstration','1'}
        name = 'Demonstration, 1';
        
        p(1).dg = 50;
        p(1).sg = 1.4;
        p(1).rhog = 12000; % density of gold-ish
        p(1).smd = 1.3;
        p(1).Dm = 3;
        type{1} = 'logn';

        p(2).dg = 200;
        p(2).sg = 1.4;
        p(2).rhog = 500; % density of salt
        p(2).smd = 1.3;
        p(2).Dm = 2.3;
        type{2} = 'logn';

    case {'soot-surrogate','2'}
        name = 'Soot surrogate, 2';
        
        p.dg = 127;
        p.sg = 1.72;
        p.rhog = 626;
        p.smd = 1.46;
        p.Dm = 2.34;
        type{1} = 'logn';

    case {'Buckley','Hogan','3'}
        name = 'Buckley-Hogan, 3';
        
        p(1).dg = 200;
        p(1).sg = 1.5;
        p(1).rhog = 10000;
        p(1).smd = 0.15;
        p(1).Dm = 3;
        type{1} = 'cond-norm';

        p(2).dg = 300;
        p(2).sg = 2.2;
        p(2).rhog = 1000;
        p(2).smd = 0.15;
        p(2).Dm = 3;
        type{2} = 'cond-norm';

    case {'narrow','4'}
        name = 'Narrow, 4';
        
        p.dg = 125;
        p.sg = 1.5;
        p.rhog = 1000;
        p.smd = 1.05;
        p.Dm = 3;
        type{1} = 'logn';
        
    otherwise
        p = [];
        type = [];
end

sm_fun = @(ii) 10^sqrt(log10(p(ii).smd)^2+...
    p(ii).Dm^2*log10(p(ii).sg)^2);
for ii=1:length(p)
    p(ii).sm = sm_fun(ii);
end

end

