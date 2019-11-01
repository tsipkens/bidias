
% GRID     Responsible for discretizing space as a grid and related operations.
% Author:  Timothy Sipkens, 2019-02-03
%-------------------------------------------------------------------------%
%   The grid class is currently used when a simple discretization of
%   two-dimensional space is required. It then takes either the span 
%   of spcae to be covered or pre-defined edge vectors to form a grid.
%   
%   See constructor method for list of other variables required 
%   for creation.
%=========================================================================%

%-- Class definition -----------------------------------------------------%
classdef Grid
    
    %-- Grid properties --------------------------------------------------%
    properties
        discrete = 'logarithmic'; % discretization to be applied to the edges
        dim = 2; % number of dimensions for mesh
        span = []; % span of values in each dimension
        
        elements = []; % contains position element centers
        nodes = []; % contains position of nodes surrounding elements
        edges = []; % vector containing edge points of element centers
        
        ne = []; % number of elements in each dimenion
        Ne = []; % total number of elements, i.e. prod(ne)
        
        nn = []; % number of nodes in each dimension
        Nn = []; % total number of nodes, i.e. prod(nn)
    end
    
    
    %-- Grid methods -----------------------------------------------------%
    methods
        %== GRID =========================================================%
        %   Class constructor.
        %-----------------------------------------------------------------%
        % Inputs: 
        %   span_edges  Either (i) a span over which discretization occurs 
        %               or (ii) a cell of edge vectors
        %   nn          If a span is specified, this is the number of nodes
        %               in each dimension
        %   discrete    Specifies type of discretization, used for
        %               marginalization and/or discretization
        %               Possible values: 'linear' or 'logarithmic' (default)
        %-----------------------------------------------------------------%
        function obj = Grid(span_edges,ne,discrete)
        
            if isa(span_edges,'cell') % consider case where edges are given
                obj.edges = span_edges;
                obj.ne = [length(span_edges{1}),...
                    length(span_edges{2})];
                obj.span = [min(span_edges{1}),max(span_edges{1});...
                    min(span_edges{2}),max(span_edges{2})];
            else % otherwise, consider case where span is given
                obj.span = span_edges;
                obj.ne = ne;
            end
            obj.Ne = prod(obj.ne);
            
            if exist('discrete','var') % if discretization scheme is specified
                if ~isempty(discrete)
                    obj.discrete = discrete;
                end
            end
            
            obj = obj.mesh; % generates grid points
        end
        %=================================================================%
        
        
        %== MESH =========================================================%
        %   Responsible for generating a mesh represented by a series of nodes.
        %   Author:	Timothy Sipkens, 2019-02-03
        %-----------------------------------------------------------------%
        %   Currently setup to do simple linear or logarithmic spaced 
        %   quadrilateral mesh.
        %
        %   obj.nodes contains the position of each of the nodes as
        %   a matrix, with a row for each node and a column for each 
        %   dimension.
        %-----------------------------------------------------------------%
        function obj = mesh(obj)
            
            %-- If required, generate edge discretization vectors --------%
            if isempty(obj.edges)
                for ii=1:obj.dim % loop through both dimensions
                    if strcmp('linear',obj.discrete)
                        obj.edges{ii} = linspace(obj.span(ii,1),obj.span(ii,2),obj.ne(ii));
                    elseif strcmp('logarithmic',obj.discrete)
                        obj.edges{ii} = logspace(...
                            log10(obj.span(ii,1)),log10(obj.span(ii,2)),obj.ne(ii));
                    end
                end
                obj.Ne = prod(obj.ne);
            end
            
            %-- Generate elements ----------------------------------------%
            [grid{1},grid{2}] = ndgrid(obj.edges{1},obj.edges{2});
            obj.elements(:,1) = grid{1}(:); % vectorize output
            obj.elements(:,2) = grid{2}(:);
            
            %-- Generate nodes -------------------------------------------%
            for ii=1:obj.dim
                if strcmp(obj.discrete,'logarithmic')
                    r_m = exp((log(obj.edges{ii}(2:end))+log(obj.edges{ii}(1:(end-1))))./2); % mean of edges
                    obj.nodes{ii} = [exp(2*log(obj.edges{ii}(1))-log(r_m(1))),...
                        r_m, exp(2*log(obj.edges{ii}(end))-log(r_m(end)))];
                elseif strcmp(obj.discrete,'linear')
                    r_m = (obj.edges{ii}(2:end)+obj.edges{ii}(1:(end-1)))./2; % mean of edges
                    obj.nodes{ii} = [2*obj.edges{ii}(1)-r_m(1),...
                        r_m, 2*obj.edges{ii}(end)-r_m(end)];
                end
            end

        end
        %=================================================================%
        
        
        %== PROJECT =========================================================%
        %   Project x onto current grid. Uses simple linear
        %   interpolation for this purpose. The parameter 'edges_old' 
        %   corresponds to edges of the original grid.
        function x = project(obj,edges_orig,x)
            
            n1 = length(edges_orig{1});
            n2 = length(edges_orig{2});
            x = reshape(x,[n1,n2]);
            
            [edges1,edges2] = ndgrid(obj.edges{1},obj.edges{2});
            [edges_old1,edges_old2] = ndgrid(edges_orig{1},edges_orig{2});
            F = griddedInterpolant(edges_old1,edges_old2,x,'linear','linear');
            
            x = F(edges1,edges2);
            x = x(:);
            
        end
        %=================================================================%
        
        
        %== DR ===========================================================%
        %   Calculates the differential area of the elements in the grid.
        function [dr,dr1,dr2] = dr(obj)
            
            dr_0 = cell(obj.dim,1);
            for ii=1:obj.dim
                if strcmp(obj.discrete,'logarithmic')
                    dr_0{ii} = log10(obj.nodes{ii}(2:end))-...
                        log10(obj.nodes{ii}(1:(end-1)));
                        
                elseif strcmp(obj.discrete,'linear')
                    dr_0{ii} = obj.nodes{ii}(2:end)-...
                        obj.nodes{ii}(1:(end-1));
                end
            end
            
            [dr1,dr2] = ndgrid(dr_0{1},dr_0{2});
            dr = dr2(:).*dr1(:);
            
        end
        %=================================================================%
        
        
        %== GRAD =========================================================%
        %   Calculates the gradient in x. Uses simple first order 
        %   differences.
        function out = grad(obj,x)
            
            [~,dr1,dr2] = obj.dr;
            
            [grad2,grad1] = gradient(reshape(x,obj.nn));
            grad1 = grad1.*dr1;
            grad2 = grad2.*dr2;
            
            out = [grad1(:),grad2(:)];
            
        end
        %=================================================================%
        
        
        %== MARGINALIZE ==================================================%
        %   Marginalized for distribution in each dimension.
        %   Uses Euler's method to integrate over domain.
        function [marg,tot] = marginalize(obj,x)
            
            [dr,dr1,dr2] = obj.dr; % generate differential area of elements
            
            tot = sum(x.*dr); % integrated total
            
            x = obj.reshape(x);
            
            marg{1} = sum(dr2.*x,2); % integrate over diameter
            marg{2} = sum(dr1.*x,1); % integrate over mass
            
        end
        %=================================================================%
        
        
        %== REBASE =======================================================%
        %   Function to evaluate uniform basis functions. Outputs matrix 
        %   to be multiplied by original A, that is to transform the kernel.
        function B = rebase(obj,obj_old)
            
            for ii=1:obj.dim
                dr_inv = 1./(obj_old.nodes{ii}(2:end)-obj_old.nodes{ii}(1:(end-1)));

                t0 = min(1,max(0,...
                    bsxfun(@times,...
                    obj_old.nodes{ii}(2:end)-obj.nodes{ii}(1:(end-1))',dr_inv)...
                    ));
                t1 = min(1,max(0,...
                    bsxfun(@times,...
                    obj.nodes{ii}(2:end)'-obj_old.nodes{ii}(1:(end-1)),dr_inv)...
                    ));
                t2{ii} = sparse(min(t0,t1));
            end
            
            ind_old_tot = 1:obj_old.Ne;
            ind1_old = mod(ind_old_tot-1,obj_old.ne(1))+1;
            ind2_old = ceil(ind_old_tot./obj_old.ne(1));
            
            ind_tot = 1:obj.Ne;
            ind1 = mod(ind_tot-1,obj.ne(1))+1;
            ind2 = ceil(ind_tot./obj.ne(1));
            
            B = (t2{1}(ind1,ind1_old).*t2{2}(ind2,ind2_old))';
            
        end
        %=================================================================%
        
        
        %== RESHAPE ======================================================%
        %   A simple function to reshape a vector based on the grid.
        function x = reshape(obj,x)
            x = reshape(x,obj.ne);
        end
        %=================================================================%
        
        
        %== VECTORIZE ====================================================%
        %   A simple function to vectorize 2D gridded data.
        %-----------------------------------------------------------------%
        % Outputs: 
        %   x	Vectorized data
        %   t1  Vectorized element centers for first dimension
        %   t2  Vectorized element centers for second dimension
        %-----------------------------------------------------------------%
        function [x,t1,t2] = vectorize(obj,x)
            if exist('x','var'); x = x(:); else; x = []; end
            
            if nargout>1; t1 = repmat(obj.edges{1}',[obj.ne(2),1]); end
            if nargout>2
                t2 = repmat(obj.edges{2},[obj.ne(1),1]);
                t2 = t2(:);
            end
        end
        %=================================================================%
        
        
        %== RAY_SUM ======================================================%
        %   Perfrom a ray sum for a given ray and the current grid. 
        %   Based on:	Code from Samuel Grauer
        %   Author:     Timothy Sipkens, 2019-07-14
        %-----------------------------------------------------------------%
        % Inputs:
        %   v0      A single point on the line
        %	slope   Slope of the line
        %   f_bar   Flag for progress bar
        % Outputs:
        %   A       Ray-sum matrix
        %-----------------------------------------------------------------%
        function A = ray_sum(obj,v0,slope,opt_bar)
            
            %-- Preallocate arrays ---------------------------------------%
            m = size(slope,1);
            A = spalloc(m,obj.Ne,0.1*m*obj.Ne); % assume 10% full
            
            %-- Compute ray-sums -----------------------------------------%
            if opt_bar, tools.textbar(); end
            for ii=1:m
                
                %-- Ray vector --%
                dv = [1,slope]; % convert slope to step vector along line
                dv = dv/norm(dv);
                dv(dv == 0) = 1e-10; % for stability
                v0(v0 == 0) = 1e-10; % for stability
                
                
                %-- Line intersections --%
                %   Use parametric representation of the line and find two
                %   intersections or each element at a time
                [~,drx,dry] = obj.dr;
                dr = [drx(:),dry(:)];
                ttmp = (log10(obj.elements(:,[2,1]))-...
                    dr./2-v0)./dv; % minimum of element
                tmax = (log10(obj.elements(:,[2,1]))+...
                    dr./2-v0)./dv; % maximum of element
                
                
                %-- Corrections --%
                %   Decide which points would correspond to transecting the
                %   pixel
                tmin = max(min(ttmp,tmax),[],2);
                tmax = min(max(ttmp,tmax),[],2);
                tmax(tmax<tmin) = tmin(tmax<tmin); % check if corsses pixel
                
                
                %-- Convert back to [x,y] --%
                rmin  = v0+tmin.*dv; % location of intersect with min. of pixel
                rmax = v0+tmax.*dv; % location of intersect with max. of pixel
                chord = sqrt(sum((rmax-rmin).^2,2)); % chord length
                
                
                %-- Ray-sum matrix --%
                [~,jj,a] = find(chord');
                if ~isempty(a)
                    A(ii,:) = sparse(1,jj,a,1,obj.Ne,0.1*obj.Ne);
                end
                if opt_bar, tools.textbar(ii/m); end
                
            end
            
        end
        %=================================================================%
        
        
        %== FIT_MASS_MOB =================================================%
        %   Fits the mass-mobility relation to the returned distribution.
        %   Author: Timothy Sipkens, 2019-07-15
        function [Dm,k,rho_100] = fit_mass_mob(obj,x,slope,opt_plot)
            
            if ~exist('slope','var'); slope = []; end
            if ~exist('opt_plot','var'); opt_plot = []; end
            
            if isempty(slope); slope = 3; end
            if isempty(opt_plot); opt_plot = 1; end
            
            %-- Get location of maximum pixel ----------------------------%
            t0 = obj.reshape(x);
            [t1,t2] = find(t0==max(max(t0)));
            v0 = log10([obj.edges{2}(t2),obj.edges{1}(t1)]);
            
            %-- Proceed with fitting -------------------------------------%
            y0 = [v0(2),slope];
            B_fun = @(y) -obj.ray_sum([v0(1),y(1)],y(2),0)*x;
            
            y1 = fminsearch(B_fun,y0);
            Dm = y1(2);
            k = 10.^(y1(1)-Dm*v0(1));
            
            rho_100 = 6/pi*k*(100^(Dm-3))*1e9;
                % expected effective density at dm = 100 nm
                % 1e9 converts from fg/nm^3 to kg/m^3
            
            if opt_plot; obj.overlay_line([v0(1),y1(1)],y1(2),'w'); end
            
        end
        %=================================================================%
        
        
        %== PLOT2D =======================================================%
        %   Plots x as a 2D function on the grid.
        %   Author: Timothy Sipkens, 2018-11-21
        function [h,x] = plot2d(obj,x)
            
            x = obj.reshape(x);
            
            if strcmp('linear',obj.discrete)
                imagesc(obj.edges{2},obj.edges{1},x);
                set(gca,'YDir','normal');
                
                xlim(obj.span(2,:));
                ylim(obj.span(1,:));
                
            elseif strcmp('logarithmic',obj.discrete)
                imagesc(log10(obj.edges{2}),log10(obj.edges{1}),x);
                set(gca,'YDir','normal');
                
                xlim(log10(obj.span(2,:)));
                ylim(log10(obj.span(1,:)));
            end

            if nargout>0; h = gca; end
            
        end
        %=================================================================%
        
        
        %== PLOT2D_MARG ==================================================%
        %   Plots x as a 2D function on the grid, with marginalized distributions.
        %   Author: Timothy Sipkens, 2018-11-21
        function [h,x_m] = plot2d_marg(obj,x,obj_t,x_t)
        
            subplot(4,4,[5,15]);
            obj.plot2d(x);
            
            x_m = obj.marginalize(x);
            
            %-- Plot marginal distribution (dim 2) -----------------------%
            subplot(4,4,[1,3]);
            marg_dim = 2;
            stairs(obj.nodes{marg_dim},...
                [x_m{marg_dim},0],'k');
            xlim([min(obj.edges{marg_dim}),max(obj.edges{marg_dim})]);
            set(gca,'XScale','log');
            
            if nargin>2 % also plot marginal of the true distribution
                x_m_t = obj_t.marginalize(x_t);
                
                hold on;
                plot(obj_t.nodes{marg_dim},...
                    [x_m_t{marg_dim},0],'color',[0.6,0.6,0.6]);
                hold off;
            end
            
            %-- Plot marginal distribution (dim 1) -----------------------%
            subplot(4,4,[8,16]);
            marg_dim = 1;
            stairs([0;x_m{marg_dim}],...
                obj.nodes{marg_dim},'k');
            ylim([min(obj.edges{marg_dim}),max(obj.edges{marg_dim})]);
            set(gca,'YScale','log');
            
            if nargin>2 % also plot marginal of the true distribution
                hold on;
                plot([0;x_m_t{marg_dim}],...
                    obj_t.nodes{marg_dim},'color',[0.6,0.6,0.6]);
                hold off;
            end
            
            subplot(4,4,[5,15]);
            if nargout>0; h = gca; end

        end
        %=================================================================%
        
        
        %== PLOT_MARGINAL ================================================%
        %   Plot marginal distributions
        %   Author:	Timothy Sipkens, 2019-07-17
        %   Note: 'x' can be a cell array containing multiple x vectors
        %-----------------------------------------------------------------%
        function [] = plot_marginal(obj,x,dim,x0)
            
            %-- Parse inputs ---------------------------------------------% 
            if ~iscell(x); x = {x}; end
                % if input is not cell, covert it to one
            
            if ~exist('dim','var'); dim = []; end
            if ~exist('x0','var'); x0 = []; end
            
            if isempty(dim); dim = 1; end
            if isempty(x0); x0 = x{1}; end
            
            x0_m = obj.marginalize(x0); % reference case
            
            
            %-- Plot entries in x ----------------------------------------%
            for ii=1:length(x) % plot other provided x
                x_m = obj.marginalize(x{ii});
                
                %-- Plot difference --%
                subplot(3,1,1);
                if ~isempty(findall(gca,'type','line')); hold on; end
                    % if not the first line in the plot, hold on
                semilogx(obj.edges{dim},x_m{dim}-x0_m{dim});
                hold off;
                
                %-- Plot marginal distribution --%
                subplot(3,1,2:3);
                if ~isempty(findall(gca,'type','line')); hold on; end
                    % if not the first line in the plot, hold on
                semilogx(obj.edges{dim},x_m{dim});
                hold off;
            end
            
            
            %-- Set axes limits ------------------------------------------%
            subplot(3,1,2:3);
            xlim([min(obj.edges{dim}),max(obj.edges{dim})]);
            
            subplot(3,1,1);
            xlim([min(obj.edges{dim}),max(obj.edges{dim})]);
        end
        %=================================================================%
        
        
        %== PLOT_CONDITIONAL =============================================%
        %   Plot conditional distributions
        %   Author:	Timothy Sipkens, 2019-07-17
        %   Note: 'x' can be a cell array containing multiple x vectors
        function [] = plot_conditional(obj,x,dim,ind,x0)
            
            %-- Parse inputs ---------------------------------------------% 
            if ~iscell(x); x = {x}; end
                % if input is not cell, covert it to one
            
            if ~exist('dim','var'); dim = []; end
            if ~exist('ind','var'); ind = []; end
            if ~exist('x0','var'); x0 = []; end
            
            if isempty(dim); dim = 1; end
            if isempty(ind); ind = round(obj.ne(dim)/2); end
            if isempty(x0); x0 = x{1}; end
            
            x0 = reshape(x0,obj.ne); % reference case
            if dim==1; x0_c = x0(:,ind); end
            if dim==2; x0_c = x0(ind,:); end
            
            
            %-- Plot entries in x ----------------------------------------%
            for ii=1:length(x) % plot other provided x
                x_c = reshape(x{ii},obj.ne);
                if dim==1; x_c = x_c(:,ind); end
                if dim==2; x_c = x_c(ind,:); end
                
                %-- Plot difference --%
                subplot(3,1,1);
                if ~isempty(findall(gca,'type','line')); hold on; end
                    % if not the first line in the plot, hold on
                semilogx(obj.edges{dim},x_c-x0_c);
                hold off;
                
                %-- Plot marginal distribution --%
                subplot(3,1,2:3);
                if ~isempty(findall(gca,'type','line')); hold on; end
                    % if not the first line in the plot, hold on
                semilogx(obj.edges{dim},x_c);
                hold off;
            end
            
            
            %-- Set axes limits ------------------------------------------%
            subplot(3,1,2:3);
            xlim([min(obj.edges{dim}),max(obj.edges{dim})]);
            
            subplot(3,1,1);
            xlim([min(obj.edges{dim}),max(obj.edges{dim})]);
        end
        %=================================================================%
        
        
        %== OVERLAY_LINE =================================================%
        %   Plots a line on top of the current grid
        %   Author:	Timothy Sipkens, 2019-07-15
        %-----------------------------------------------------------------%
        % Inputs:
        %   r0      A single point on the line
        %	slope   Slope of the line
        %   c_spec  Color specification string, e.g. 'k' for a black line
        % Outputs:
        %   h       Line object
        %-----------------------------------------------------------------%
        function h = overlay_line(obj,r0,slope,cspec)
            
            if ~exist('cspec','var'); cspec = 'w'; end
            
            rmin = log10(min([obj.edges{:}]));
            rmax = log10(max([obj.edges{:}]));
            
            hold on;
            h = plot([rmin,rmax],...
                [r0(2)+slope*(rmin-r0(1)),...
                r0(2)+slope*(rmax-r0(1))],cspec);
            hold off;
            
            if nargout==0; clear h; end
            
        end
        %=================================================================%
        
    end
    
    
    methods(Static)    
        %== EFF_IND ======================================================%
        %   Returns the effective index between two grids, that is the 
        %   index including fraction.
        function ind = eff_ind(edges,r)
            
            dim = length(edges); % will be number of dimensions
            
            for ii=1:dim
                ind_bool = and(...
                    r(:,ii)>=edges{ii}(1:(end-1)),...
                    r(:,ii)<edges{ii}(2:end));
                
                ind_temp = sum(cumprod(ind_bool==0,2),2)+1;
                ind_temp(r(:,ii)==edges{ii}(end)) = ...
                    ind_temp(r(:,ii)==edges{ii}(end))-1;
                        % condition for points at end of domain
                ind(:,ii) = ind_temp;
            end
        end
        %=================================================================%
    end
end

