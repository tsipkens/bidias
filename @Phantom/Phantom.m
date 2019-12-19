
% PHANTOM  Class containing properties and methods for storing phantoms.
% Author:  Timothy Sipkens, 2019-07-08
%=========================================================================%

classdef Phantom

    %-- Phantom properties -----------------------------------------------%
    properties
        type = [];  % optional name for the phantom
        modes = {}; % types of distribution for each mode
        n_modes = []; % number of modes

        mu = [];    % center of bivariate distribution
        Sigma = []; % covariance of bivariate distribution
        R = [];     % correlation matrix

        x = [];     % phantom evaluated on default grid
        grid = [];  % default grid the phantom is to be represented
                    % on generally a high resolution instance of the
                    % Grid class

        p struct = struct();
                    % parameters relevant to mass-mobility distributions
                    % fields include mg, sm, sd, rho_100, etc.
    end
    


    methods
        %== PHANTOM ======================================================%
        %   Intialize phantom object.
        %-----------------------------------------------------------------%
        % Inputs:
        %   type        The type of phantom specified as a string
        %               (e.g. 'standard', 'mass-mobility', '1')
        %   span_grid   Either (i) the span over which the phantom is to be
        %               evaluated or (ii) a grid on which the phantom is to
        %               be evaluated
        %   mu_p        Either a vector of distribution means or a vector p
        %               specifying the mass-mobility parameters
        %   Sigma_modes Either the covariance matrix for the distribution
        %               or the number of modes in the distribution
        %               (e.g. 'logn','cond-norm')
        %-----------------------------------------------------------------%
        function [obj] = Phantom(type,span_grid,mu_p,Sigma_modes)

            %-- Parse inputs ---------------------------------------------%
            if nargin==0; return; end % return empty phantom

            if ~exist('span_grid','var'); span_grid = []; end
            if isempty(span_grid); span_grid = [10^-1.5,10^1.5;20,10^3]; end
            %-------------------------------------------------------------%

            %== Assign parameter values - 3 options ======================%
            switch type

                %-- OPTION 1: Standard bivariate lognormal distribution --%
                case {'standard'}
                    n_modes = length(mu_p);
                    if ~iscell(mu_p); n_modes = 1; end

                    obj.mu = mu_p;
                    obj.Sigma = Sigma_modes;

                    [obj.modes{1:n_modes}] = 'logn';

                    p = obj.cov2p(obj.mu,obj.Sigma,obj.modes);
                        % mass-mobility equivlanet params.
                    obj.p = obj.fill_p(p); % get mg as well

                %-- OPTION 2: Using a mass-mobility parameter set (p) ----%
                case {'mass-mobility'} % for custom mass-mobility phantom
                                       % sepcified using a p structure
                    obj.type = 'mass-mobility';
                    obj.modes = Sigma_modes;
                    obj.p = mu_p;

                    if ~any(strcmp('cond-norm',Sigma_modes))
                        [obj.mu,obj.Sigma] = obj.p2cov(obj.p,obj.modes);
                    end

                %-- OPTION 3: Use a preset or sample distribution --------%
                otherwise % check if type is a preset phantom
                    [p,modes,type] = obj.preset_phantoms(type);
                    if isempty(p); error('Invalid phantom call.'); end

                    obj.type = type;
                    obj.modes = modes;

                    obj.p = obj.fill_p(p);

                    if ~any(strcmp('cond-norm',modes))
                        [obj.mu,obj.Sigma] = obj.p2cov(obj.p,obj.modes);
                    end
            end

            obj.n_modes = length(obj.modes); % get number of modes


            %-- Generate a grid to evaluate phantom on -------------------%
            if isa(span_grid,'Grid') % if grid is specified
                obj.grid = span_grid;
            else % if span is specified, create grid
                n_t = [540,550]; % resolution of phantom distribution
                obj.grid = Grid(span_grid,...
                    n_t,'logarithmic'); % generate grid of which to represent phantom
            end


            %-- Evaluate phantom -----------------------------------------%
            if any(strcmp('cond-norm',obj.modes))
                obj.x = obj.eval_p(obj.p);
                    % special evaluation for conditional normal conditions
            else
                obj.x = obj.eval(obj.mu,obj.Sigma);
            end
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
        %   Plots the phantom distribution.
        %   Author:     Timothy Sipkens, 2019-07-08
        function [h] = plot(obj)
            h = obj.grid.plot2d_marg(obj.x);
        end
        %=================================================================%



        %== EVAL =========================================================%
        %   Generates a distribution from the phantom mean and covariance.
        %   Author:  Timothy Sipkens, 2019-10-29
        %   NOTE: Does not work for conditionally-normal distributions
        %         (which cannot be defined with mu and Sigma).
        function [x] = eval(obj,mu,Sigma)

            if ~iscell(mu); mu = {mu}; end
            if ~iscell(Sigma); Sigma = {Sigma}; end

            [~,m_vec,d_vec] = obj.grid.vectorize();

            %-- Assign other parameters of distribution ------------------%
            x = zeros(size(m_vec));
            for ll=1:obj.n_modes % loop through distribution modes
                x = x + mvnpdf(log10([m_vec,d_vec]),mu{ll},Sigma{ll});
            end

            %-- Reweight modes -------------------------------------------%
            x = x./obj.n_modes;
        end
        %=================================================================%



        %== EVAL_P =======================================================%
        %   Generates a distribution from p as required for conditionally-
        %   normal modes.
        %   Author:  Timothy Sipkens, 2019-10-29
        function [x] = eval_p(obj,p)

            [~,m_vec,d_vec] = obj.grid.vectorize();

            m_fun = @(d,ll) log(p(ll).m_100.*((d./100).^p(ll).Dm));
                % geometric mean mass in fg as a function of d

            %-- Evaluate phantom mass-mobility distribution ---------------%
            x = zeros(length(m_vec),1); % initialize distribution parameter
            for ll=1:obj.n_modes % loop over different modes
                if strcmp(obj.modes{ll},'logn')
                    p_m = lognpdf(m_vec,m_fun(d_vec,ll),log(p(ll).smd));
                else
                    p_m = normpdf(m_vec,...
                        exp(m_fun(d_vec,ll)),p(ll).smd.*...
                        exp(m_fun(d_vec,ll)));
                end

                p_temp = lognpdf(d_vec,log(p(ll).dg),log(p(ll).sg)).*p_m;
                x = x+p_temp;
            end

            %-- Reweight modes and transform to log-log space ------------%
            x = x./obj.n_modes;
            x = x.*(d_vec.*m_vec).*log(10).^2;
                % convert to [log10(m),log10(d)]T space
        end
        %=================================================================%



        %== MASS2RHO =====================================================%
        %   Convert a mass-mobility phantom to an effective density-mobility
        %   phanatom. Output is a new phantom in the transformed space.
        %   Author:  Timothy Sipkens, 2019-10-31
        function [phantom] = mass2rho(obj,grid_rho)

            A = [1,-3;0,1]; % corresponds to mass-mobility relation

            mu_rhod = (A*obj.mu'+[log10(6/pi)+9;0])';
            Sigma_rhod = A*obj.Sigma*A';

            phantom = Phantom('standard',grid_rho,mu_rhod,Sigma_rhod);

        end
        %=================================================================%
    end



    methods (Static)
        %== PRESET_PHANTOMS ==============================================%
        % Returns a set of parameters for preset/sample phantoms.
        % Definition is in an external function.
        [p,modes,type] = preset_phantoms(obj,type);
        %=================================================================%


        %== FILL_P =======================================================%
        %   Generates the remainder of the components of p.
        %   Author:  Timothy Sipkens, 2019-10-30
        function p = fill_p(p)

            n_modes = length(p);

            %-- Assign other parameters of distribution ------------------%
            for ll=1:n_modes % loop through distribution modes
                p(ll).mg = 1e-9*p(ll).rhog*pi/6*...
                    (p(ll).dg^3); % geometric mean mass in fg

                if ~isfield(p,'rhog'); p.rho = []; end

                if ~isempty(p(ll).rhog) % use effective density at dg
                    p(ll).rho_100 = p(ll).rhog*(100/p(ll).dg)^(p(ll).Dm-3);
                end

                p(ll).m_100 = 1e-9*p(ll).rho_100*pi/6*100^3;
                p(ll).rhog = p(ll).rho_100*((p(ll).dg/100)^(p(ll).Dm-3));
                p(ll).k = p(ll).m_100/(100^p(ll).Dm);
            end

        end
        %=================================================================%



        %== VEC2P ========================================================%
        %   Function to format phantom parameters from a vector, t.
        %   Author:  Timothy Sipkens, 2019-07-18
        function [p] = vec2p(vec,n_modes)
            t_length = length(vec);
            n = n_modes*t_length;

            p.dg = vec(1:t_length:n);
            p.sg = vec(2:t_length:n);
            p.rho_100 = vec(3:t_length:n);
            p.smd = vec(4:t_length:n);
            p.Dm = vec(5:t_length:n);
        end
        %=================================================================%



        %== COV2P ========================================================%
        %   Function to convert covariance matrix and mean to p.
        %   Author:  Timothy Sipkens, 2019-10-29
        function [p] = cov2p(mu,Sigma,modes)

            if ~iscell(mu); mu = {mu}; end
            if ~iscell(Sigma); Sigma = {Sigma}; end

            p = [];
            for ll=length(mu):-1:1 % loop through modes
                p(ll).dg = 10.^mu{ll}(2);
                p(ll).mg = mu{ll}(1);

                if strcmp(modes{ll},'logn') % if lognormal distribution, convert to geometric mean
                    p(ll).mg = 10.^p(ll).mg;
                end

                p(ll).sg = 10^sqrt(Sigma{ll}(2,2));
                p(ll).sm = 10^sqrt(Sigma{ll}(1,1));

                Sigma_inv = inv(Sigma{ll});
                p(ll).smd = 10^sqrt(1/Sigma_inv(1,1));

                p(ll).Dm = Sigma{ll}(1,2)/Sigma{ll}(2,2);
                    % corresponds to slope of "locus of vertical"
                    % (Friendly, Monette, and Fox, 2013)
                    % can be calculated as Dm = corr*sy/sx

                % t0 = eigs(rot90(Sigma{ll},2),1);
                % p(ll).ma = (t0-Sigma{ll}(2,2))./Sigma{ll}(1,2);
                    % calculate the major axis slope

                p(ll).rhog = p(ll).mg/(pi*p(ll).dg^3/6)*1e9;
            end

            p = Phantom.fill_p(p);
        end
        %=================================================================%



        %== P2COV ========================================================%
        %   Function to convert p to a covariance matrix and mean.
        %   Author:  Timothy Sipkens, 2019-10-29
        function [mu,Sigma] = p2cov(p,modes)

            for ll=length(p):-1:1
                mu{ll} = log10([p(ll).mg,p(ll).dg]);
                    % use geometric mean

                if strcmp(modes{ll},'logn')
                    Sigma{ll} = inv([(1/log10(p(ll).smd))^2,...
                        -p(ll).Dm/log10(p(ll).smd)^2;...
                        -p(ll).Dm/log10(p(ll).smd)^2,...
                        1/log10(p(ll).sg)^2+p(ll).Dm^2/log10(p(ll).smd)^2]);
                else
                    Sigma{ll} = inv([(1/p(ll).smd)^2,...
                        -p(ll).Dm/p(ll).smd^2;...
                        -p(ll).Dm/p(ll).smd^2,...
                        1/log10(p(ll).sg)^2+p(ll).Dm^2/p(ll).smd^2]);
                end
            end

            if length(mu)==1
                mu = mu{1};
                Sigma = Sigma{1};
            end
        end
        %=================================================================%



        %== FIT ==========================================================%
        %   Fits a phantom to a given set of data, x, defined on a given
        %   grid. Outputs a fit phantom object.
        %-----------------------------------------------------------------%
        % Inputs:
        %   x - input data (2D distribution data)
        %   grid - 'Grid' object on which input data is evaluated
        %-----------------------------------------------------------------%
        function phantom = fit(x,grid)

            [~,vec1,vec2] = grid.vectorize();

            corr2cov = @(sigma,R) diag(sigma)*R*diag(sigma);

            fun_pha = @(y) y(1).*mvnpdf(log10([vec1,vec2]),[y(2),y(3)],...
                corr2cov([y(4),y(5)],[1,y(6);y(6),1]));
            y0 = [max(x),0,2.3,0.3,0.2,0.99];
                % [C,mg,dg,sm,sg,corr]

            y1 = lsqnonlin(@(y) fun_pha(y)-x, y0, ...
                [0,-10,-10,0,0,-1],[inf,10,10,10,3,1]);

            mu = [y1(2),y1(3)];
            sigma = [y1(4),y1(5)];
            Sigma = corr2cov(sigma,[1,y1(6);y1(6),1]);

            phantom = Phantom('standard',grid,mu,Sigma);
            phantom.type = 'standard-fit';
        end
        %=================================================================%
    end

end
