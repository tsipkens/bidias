
classdef Grid
% GRID Responsible for meshing and marginalization.
% Author:   Timothy Sipkens, 2019-02-03
%   
%-------------------------------------------------------------------------%
%   The grid class is currently used when a simple discretization of
%   space of two-dimensional space is required. 
%   
%   See constructor method for list of variables required for creation.
%   
%   This class currently only generates rectangular meshes with
%   uniform or logarithmic spacing or specified edge vectors.
%-------------------------------------------------------------------------%
    
    
    properties
        discrete = 'logarithmic'; % discretization to be applied to the edges
        dim = []; % number of dimensions for mesh
        span = []; % span of values in each dimension
        
        elements = []; % contains position element centers
        edges = []; % vector containing edge points of element centers
        nodes = []; % contains position of nodes surrounding elements
        
        ne = []; % number of elements in each dimenion
        Ne = []; % total number of elements
        
        nn = []; % number of nodes in each dimension
        Nn = []; % total number of nodes
    end
    
    
    methods
        function obj = Grid(span_edges,ne,discrete)
        % GRID Constructor method for this class.
        % Author:   Timothy Sipkens, 2019-02-03
        %
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
        
            
            if isa(span_edges,'cell') % consider case where edges are given
                obj.edges = span_edges;
                obj.ne = [length(span_edges{1}),...
                    length(span_edges{2})];
                obj.span = [min(span_edges{1}),max(span_edges{1});...
                    min(span_edges{2}),max(span_edges{2})];
                obj.dim = 2;
                
            else % consider case where span is given
                obj.span = span_edges;
                obj.ne = ne;
                obj.dim = 2;
                
            end
            
            if exist('discrete','var') % if discretization scheme is specified
                if ~isempty(discrete)
                    obj.discrete = discrete;
                end
            end
            
            obj = obj.mesh; % generates grid points
        end
        
        
        function obj = mesh(obj)
        % MESH Responsible for generating a mesh represented by a series of nodes.
        % Author:   Timothy Sipkens, 2019-02-03
        %
        %-----------------------------------------------------------------%
        %   Currently setup to do simple linear or logarithmic spaced 
        %   quadrilateral mesh.
        %
        %   obj.nodes contains the position of each of the nodes as
        %   a matrix, with a row for each node and a column for each 
        %   dimension.
        %-----------------------------------------------------------------%
            
            
            obj.Ne = prod(obj.ne);
            
            %-- If required, generate edge discretization vectors --------%
            if isempty(obj.edges)
                for ii=1:obj.dim
                    if strcmp('linear',obj.discrete)
                        obj.edges{ii} = linspace(obj.span(ii,1),obj.span(ii,2),obj.ne(ii));
                    elseif strcmp('logarithmic',obj.discrete)
                        obj.edges{ii} = logspace(log10(obj.span(ii,1)),log10(obj.span(ii,2)),obj.ne(ii));
                    end
                end
                obj.Ne = prod(obj.ne);
            end
            
            
            %-- Generate elements ----------------------------------------%
            [grid{1},grid{2}] = ndgrid(obj.edges{1},obj.edges{2});
            obj.elements(:,1) = grid{1}(:); % vectorize output
            obj.elements(:,2) = grid{2}(:);
            
            
            %-- Generate nodes -------------------------------------------%
            obj.nodes = cell(obj.dim,1);
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
        
        
        function [h,x] = plot2d(obj,x)
        % PLOT2D Plots x on the grid, with output depending on discretization scheme.
        % Author:   Timothy Sipkens, 2018-11-21
            
            x = reshape(x,obj.ne);

            if strcmp('linear',obj.discrete)
                h = imagesc(obj.edges{2},obj.edges{1},x);
                % set(h, 'EdgeColor', 'none');
                set(gca,'YDir','normal');
            elseif strcmp('logarithmic',obj.discrete)
                h = imagesc(log10(obj.edges{2}),log10(obj.edges{1}),x);
                % set(h, 'EdgeColor', 'none');
                set(gca,'YDir','normal');
            end
        end
        
        
        function x = project(obj,edges,x)
        % PROJECT Project x onto current grid.
        % Uses simple linear interpolation for this purpose.
            
            n1 = length(edges{1});
            n2 = length(edges{2});
            x = reshape(x,[n1,n2]);
            
            [edges1,edges2] = ndgrid(obj.edges{1},obj.edges{2});
            [edges_old1,edges_old2] = ndgrid(edges{1},edges{2});
            F = griddedInterpolant(edges_old1,edges_old2,x,'linear','linear');
            
            x = F(edges1,edges2);
            x = x(:);
            
        end
        
        
        function [dr,dr1,dr2] = dr(obj)
        % DR Calculates the differential area of the elements in the grid.
            
            dr_0 = cell(obj.dim,1);
            for ii=1:obj.dim
                if strcmp(obj.discrete,'logarithmic')
                    dr_0{ii} = log(obj.nodes{ii}(2:end))-log(obj.nodes{ii}(1:(end-1)));
                elseif strcmp(obj.discrete,'linear')
                    dr_0{ii} = obj.nodes{ii}(2:end)-obj.nodes{ii}(1:(end-1));
                end
            end
            
            [dr1,dr2] = ndgrid(dr_0{1},dr_0{2});
            dr = dr2(:).*dr1(:);
            
        end
        
        
        function out = grad(obj,x)
        % GRAD Calculates the gradient in x.
        % Uses simple first order differences.
            
            [~,dr1,dr2] = obj.dr;
            
            [grad2,grad1] = gradient(reshape(x,obj.nn));
            grad1 = grad1.*dr1;
            grad2 = grad2.*dr2;
            
            out = [grad1(:),grad2(:)];
            
        end
        
        
        function [marg,tot] = marginalize(obj,x)
        % MARGINALIZE Marginalized for distribution in each dimension.
        % Uses Euler's method to integrate over domain.
            
            [dr,dr1,dr2] = obj.dr;
            
            tot = sum(x.*dr); % integrated total
            
            x = reshape(x,obj.ne);
            
            marg{1} = sum(dr2.*x,2); % integrate over diameter
            marg{2} = sum(dr1.*x,1); % integrate over mass
            
        end
        
        
        function B = phi(obj,obj_old)
        % PHI Function to evaluate uniform basis functions.
            
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
        
    end
    
    
    methods(Static)        
        function ind = eff_ind(edges,r)
            % EFF_IND Returns the effective index between two grids, that
            % is the index including fraction.
            
            dim = length(edges); % will be number of dimensions
            ind_start = zeros(length(r(:,1)),dim);
            ind_end = zeros(length(r(:,1)),dim);
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
    end
end

