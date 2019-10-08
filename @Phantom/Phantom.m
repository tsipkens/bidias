
% PHANTOM  Class containing properties and methods for storing phantoms.
% Author:  Timothy Sipkens, 2019-07-08
%=========================================================================%

%-- Class definition -----------------------------------------------------%
classdef Phantom
    
    %-- Phantom properties -----------------------------------------------%
    properties
        name = []; % optional name for the phantom
        
        p struct = struct();
        %{
            ... % struct defining distribution parameterization
            'dg',[],... % mean mobility diameter [nm]
            'sg',[],... % standard deviation of mobility diameter
            'rho',[],... % effective density at mean mobility diameter [kg/m3]
            'rho_100',[],... % effective density at d = 100 nm [kg/m3]
            'sm',[],... % standard deviation of mas
            'Dm',[],... % mass-mobility exponent
            'mg',[],...  % conditional mode of the mass distribution
            'k',[],... % mass-mobility pre-exponential factor
            'C',[],... % factor used to scale fits to the data (empty otherwise)
            'type_m',{}); % form of probability distribution for mass ('logn' or 'norm')
        %}
        
        n_modes = []; % number of modes
        % mg_fun = []; % ridge of the mass-mobility distribution (function handle)
        % rho_fun = []; % ridge of effective density-mobility distribution (function handle)
        
        x = []; % evaluated phantom
        grid = []; % grid phantom is represented on
                   % generally a high resolution mesh
    end
    
    
    %-- Phantom methods --------------------------------------------------%
    methods
        %== PHANTOM ======================================================%
        %   Intialize phantom object.
        function [obj,x] = Phantom(name,span_grid,p)
            
            %-- Assign parameter values ----------------------------------%
            switch name
                case {'demonstration','1'}
                    obj.name = 'demonstration';
                    
                    p(1).dg = 50;
                    p(1).sg = 1.4;
                    p(1).rho = 12000; % density of gold-ish
                    p(1).sm = 1.3;
                    p(1).Dm = 3;
                    p(1).type_m = 'logn';
                    
                    p(2).dg = 200;
                    p(2).sg = 1.4;
                    p(2).rho = 500; % density of salt
                    p(2).sm = 1.3;
                    p(2).Dm = 2.3;
                    p(2).type_m = 'logn';
                
                case {'2-old'}
                    obj.name = 'old-soot-surrogate';
                    
                    p.dg = 125;
                    p.sg = 1.6;
                    p.Dm = 2.3;
                    p.sm = 1.5;
                    p.k = 9400;
                    p.rho = 6*p.k/...
                        pi*p.dg^(p.Dm-3);
                    p.type_m = 'logn';
                    
                case {'soot-surrogate','2'}
                    obj.name = 'soot-surrogate';
                    
                    p.dg = 127;
                    p.sg = 1.72;
                    p.rho = 626;
                    p.sm = 1.46;
                    p.Dm = 2.34;
                    p.type_m = 'logn';
                    
                case {'Buckley','Hogan','3'}
                    obj.name = 'Buckley';
                    
                    p(1).dg = 200;
                    p(1).sg = 1.5;
                    p(1).rho = 10000;
                    p(1).sm = 0.15;
                    p(1).Dm = 3;
                    p(1).type_m = 'norm';
                    
                    p(2).dg = 300;
                    p(2).sg = 2.2;
                    p(2).rho = 1000;
                    p(2).sm = 0.15;
                    p(2).Dm = 3;
                    p(2).type_m = 'norm';
                    
                case {'narrow','4'}
                    obj.name = 'narrow';
                    
                    p.dg = 125;
                    p.sg = 1.5;
                    p.rho = 1000;
                    p.sm = 1.05;
                    p.Dm = 3;
                    p.type_m = 'logn';
                    
                otherwise % for custom phantom
                    if ~exist('p','var'); error('Specify phantom.'); end
                    if isempty(p); error('Specify phantom.'); end
                    
                    obj.name = 'custom';
            end
            
            
            %-- Evaluate additional parameters ---------------------------%
            obj.n_modes = length(p); % get number of modes
            obj.p = p; % assign distribution parameters to the instance of this class
            
            
            %-- Generate a grid to evaluate phantom on -------------------%
            if isa(span_grid,'Grid') % if grid is specified
                obj.grid = span_grid;
                
            else % if span is specified
                n_t = [540,550]; % resolution of phantom distribution
                obj.grid = Grid(span_grid,... 
                    n_t,'logarithmic'); % generate grid of which to represent phantom
            end
            
            
            %-- Evaluate phantom -----------------------------------------%
            [obj.x,obj.p] = ...
                Phantom.eval_phantom_mm(p,obj.grid);
                    % evaluate phantom, updating x and p
        end
        %=================================================================%
        
        
        %== MG_FUN =======================================================%
        %   Function to evaluate mg as a function of mobility diameter.
        %   Calculation is based on the empirical mass-mobility relation.
        %   Author:  Timothy Sipkens, 2019-07-19
        function mg = mg_fun(obj,d)
            
            d_size = size(d);
            if d_size(2)>d_size(1); d = d'; d_size = size(d); end
                % transpose the mobility if necessary
            
            mg = zeros(d_size(1),obj.n_modes);
            rho = obj.rho_fun(d);
            for ll=1:obj.n_modes
                mg(:,ll) = 1e-9.*rho(:,ll).*pi./6.*(d.^3);
                    % output in fg
            end
            
        end
        %=================================================================%
        
        
        %== RHO_FUN ======================================================%
        %   Function to evaluate the effective density as a function of
        %   mobility diameter.
        %   Author:  Timothy Sipkens, 2019-07-19
        function rho = rho_fun(obj,d)
            
            rho = zeros(length(d),obj.n_modes);
            for ll=1:obj.n_modes
                rho(:,ll) = 6*obj.p(ll).k./(pi.*d.^(3-obj.p(ll).Dm));
                    % output in kg/m3
            end
            
        end
        %=================================================================%
        
        
        %== PLOT =========================================================%
        %   Plots the phantom mass-mobiltiy distribution phantom.
        %   Author:     Timothy Sipkens, 2019-07-08
        function [h] = plot(obj)
            h = obj.grid.plot2d_marg(obj.x);
        end
        %=================================================================%
        
    end
    
    methods (Static)
        %== EVAL_PHANTOM_MM ==============================================%
        %   Generates a mass-mobiltiy distribution phantom from a set 
        %   of input parameters, p.
        %   Author:  Timothy Sipkens, 2018-12-04
        function [x,p] = eval_phantom_mm(p,grid)
            
            m_vec = grid.elements(:,1);
            d_vec = grid.elements(:,2);
            n_modes = length(p);
            
            %-- Assign other parameters of distribution ------------------%
            for ll=1:n_modes % loop through distribution modes
                if ~isfield(p,'rho'); p.rho = []; end
                
                if ~isempty(p(ll).rho) % use effective density at dg
                    p(ll).rho_100 = p(ll).rho*(100/p(ll).dg)^(p(ll).Dm-3);
                end
                
                p(ll).m_100 = 1e-9*p(ll).rho_100*pi/6*100^3;
                p(ll).rho = p(ll).rho_100*((p(ll).dg/100)^(p(ll).Dm-3));
                p(ll).mg = 1e-9*p(ll).rho*pi/6*...
                    (p(ll).dg^3); % geometric mean mass in fg
                
                p(ll).k = p(ll).m_100/(100^p(ll).Dm);
            end
            
            m_fun = @(d,ll) log(p(ll).m_100.*((d./100).^p(ll).Dm));
                % geometric mean mass in fg as a function of d
            
            
            %-- Evaluate phantom mass-mobility distribution ---------------%
            x = zeros(length(m_vec),1); % initialize distribution parameter
            for ll=1:n_modes % loop over different modes
                if strcmp(p(ll).type_m,'norm')
                    p_m = normpdf(m_vec,...
                        exp(m_fun(d_vec,ll)),p(ll).sm.*...
                        exp(m_fun(d_vec,ll)));
                else
                    p_m = lognpdf(m_vec,m_fun(d_vec,ll),log(p(ll).sm));
                end
                
                p_temp = lognpdf(d_vec,log(p(ll).dg),log(p(ll).sg)).*p_m;
                x = x+p_temp;
            end
            
            
            %-- Reweight modes and transform to log-log space ------------%
            x = x./n_modes;
            x = x.*(d_vec.*m_vec).*log(10).^2;
                % convert to [log10(m),log10(d)]T space  
        end
        %=================================================================%
        
        
        %== VEC2P ========================================================%
        %   Function to format phantom parameters from a vector, t.
        %   Author:  Timothy Sipkens, 2019-07-18
        function [p] = vec2p(vec,type,n_modes)
            if ~exist('type_m','var'); type = []; end
            if isempty(type); type = 'logn'; end
            
            t_length = length(vec);
            n = n_modes*t_length;
            
            p.dg = vec(1:t_length:n);
            p.sg = vec(2:t_length:n);
            p.rho_100 = vec(3:t_length:n);
            p.sm = vec(4:t_length:n);
            p.Dm = vec(5:t_length:n);
            p.type_m = type;
        end
        %=================================================================%
        
        
        %== PHANTOM ======================================================%
        function phantom = fit(grid,x,type)
            % Inputs:
            %   grid - 'Grid' object on which input data is evaluated
            %   x - input data (2D distribution data)
            %   type - type of distr. (e.g. 'logn')
            %-----------------------------------------------------%
            
            n_modes = 1;
            
            t0 = [160,1.8,700,1.5,2.6,2.5]; % initial guess
            % t0 = [130,1.8,800,1.8,2.38,2]; % initial guess
                % format: [dg,sg,rho_100,sm,Dm,log10(C)]
            fun = @(t) (10.^t(6)).*x-...
                Phantom.eval_phantom_mm(Phantom.vec2p(t,type,n_modes),grid);
            
            t1 = lsqnonlin(fun,t0); % fit phantom to provided data
            
            p = Phantom.vec2p(t1,type,n_modes);
            
            phantom = Phantom('custom',grid,p);
            phantom.name = 'fit';
        end
        %=================================================================%
    end
    
end

