
% GRID  Responsible for discretizing space as a grid and related operations.
% Author:  Timothy Sipkens, 2019-02-03
% 
% The grid class is currently used when a simple discretization of
% two-dimensional space is required. It then takes either the span
% of spcae to be covered or pre-defined edge vectors to form a grid.
%
% See constructor method for list of other variables required
% for creation.
%=========================================================================%

classdef Grid


properties
    discrete = 'logarithmic';
                % type discretization to be applied to the edges
                % ('logarithmic' or 'linear')
    
    dim = 2;    % number of dimensions of mesh
    
    span = [];  % span covered by the grid in each dimension
                % Span applies to the center of the elements, i.e.
                % span(1,1) is the center of the first element for the
                % first dimension.
    
    ne = [];    % number of pixels/elements in each dimenion
    Ne = [];    % total number of pixels/elements, initially Ne = prod(ne)
                % reduced for partial grids
    
    edges = []; % cell of vectors containing edge points of pixel/element
                % centers, with one cell entry per dimension
    nodes = []; % contains position of nodes surrounding elements for each dimension
                % each cell has a vector of size (ne + 1).
    
    elements = [];  % contains position of pixel/element centers as a (ne x 2) vector
    nelements = []; % position of pixel/elements edges as a (ne x 4) vector
                    % [dim1_low,dim1_high,dim2_low,dim2_high]
                
    adj = [];   % adjacency matrix
    
    
    %-- Partial grid properties --------------------------------------%
    ispartial = 0; % toggle of whether the grid is 'partial' or sparse,
                   % that is having some grid points missing
                   
    missing = [];  % global indices of missing grid points for partial grids
    cut = [];      % [y-intercept,slope] used to cut partial grid
    %-----------------------------------------------------------------%
end


methods
    %== GRID =========================================================%
    %   Class constructor.
    %-----------------------------------------------------------------%
    % Inputs:
    %   span_edges  Either:
    %                (i) a span over which discretization occurs or
    %                (ii) a cell of edge vectors
    %   ne          If a span is specified, this is the number of
    %               elements/pixels in each dimension
    %   discrete    Specifies type of discretization, used for
    %               marginalization and/or discretization
    %               Possible values: 'linear' or 'logarithmic' (default)
    %-----------------------------------------------------------------%
    function obj = Grid(span_edges,ne,discrete)
        
        %-- Parse inputs ---------------------------------------------%
        if nargin==0; return; end % return empty grid
        
        if ~exist('discrete','var'); discrete = []; end
        if isempty(discrete); discrete = 'logarithmic'; end
            % if discrete is not specified, use logarithmic
        %-------------------------------------------------------------%
        
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
        obj = obj.adjacency; % get adjacency matrix
    end
    %=================================================================%
    
    

    %== MESH =========================================================%
    %   Responsible for generating a mesh represented by a series of elements.
    %   Author:	Timothy Sipkens, 2019-02-03
    %   
    %   Currently setup to do simple linear or logarithmic spaced
    %   quadrilateral mesh.
    %
    %   obj.nodes contains the position of each of the nodes as
    %   a matrix, with a row for each node and a column for each
    %   dimension.
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

        %-- Generate nodes -------------------------------------------%
        for ii=1:obj.dim
            if strcmp(obj.discrete,'logarithmic')
                r_m = exp((log(obj.edges{ii}(2:end))+...
                    log(obj.edges{ii}(1:(end-1))))./2); % mean of edges

                obj.nodes{ii} = [exp(2*log(obj.edges{ii}(1))-log(r_m(1))),...
                    r_m, exp(2*log(obj.edges{ii}(end))-log(r_m(end)))];

            elseif strcmp(obj.discrete,'linear')
                r_m = (obj.edges{ii}(2:end)+...
                    obj.edges{ii}(1:(end-1)))./2; % mean of edges

                obj.nodes{ii} = [2*obj.edges{ii}(1)-r_m(1),...
                    r_m, 2*obj.edges{ii}(end)-r_m(end)];
            end
        end

        %-- Generate vectorized lists of elements --------------------%
        %   One column per dimension
        [vec1{1},vec1{2}] = ndgrid(obj.edges{1},obj.edges{2});
        obj.elements(:,1) = vec1{1}(:); % vectorize output
        obj.elements(:,2) = vec1{2}(:);
        
        [vec1{1},vec1{2}] = ndgrid(obj.nodes{1}(1:(end-1)),...
            obj.nodes{2}(1:(end-1)));
        [vec2{1},vec2{2}] = ndgrid(obj.nodes{1}(2:end),...
            obj.nodes{2}(2:end));
        obj.nelements = [vec1{1}(:),vec2{1}(:),vec1{2}(:),vec2{2}(:)];
    end
    %=================================================================%
    
    
    
    %== ADJACENCY ====================================================%
    %   Compute the adjacency matrix for the full grid,
    %   using a four-point stencil.
    function [obj,adj] = adjacency(obj)
        
        ind1 = ones(3*prod(obj.ne),1);
        ind2 = ones(3*prod(obj.ne),1);
        vec = zeros(3*prod(obj.ne),1);
        ll = 0;
        
        for jj=1:prod(obj.ne)
            if ~(mod(jj,obj.ne(1))==0) % up pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1;
                vec(ll) = 1;
            end
            
            if ~(mod(jj-1,obj.ne(1))==0) % down pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1;
                vec(ll) = 1;
            end
            
            if jj>obj.ne(1) % left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-obj.ne(1);
                vec(ll) = 1;
            end
            
            if jj<=(prod(obj.ne)-obj.ne(1)) % right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+obj.ne(1);
                vec(ll) = 1;
            end
        end
        
        adj = sparse(ind1,ind2,vec,...
            prod(obj.ne),prod(obj.ne));
        
        obj.adj = adj;
    end
    %=================================================================%
    
    
    
    %== ADJACENCY8 ===================================================%
    %   Compute the adjacency matrix for the full grid,
    %   using an eight-point stencil.
    function [adj] = adjacency8(obj)
        
        ind1 = ones(7*prod(obj.ne),1);
        ind2 = ones(7*prod(obj.ne),1);
        vec = zeros(7*prod(obj.ne),1);
        ll = 0;
        
        for jj=1:prod(obj.ne)
            if ~(mod(jj,obj.ne(1))==0) % up pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1;
                vec(ll) = 1;
            end
            
            if ~(mod(jj-1,obj.ne(1))==0) % down pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1;
                vec(ll) = 1;
            end
            
            if jj>obj.ne(1) % left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-obj.ne(1);
                vec(ll) = 1;
            end
            
            if jj<=(obj.Ne-obj.ne(1)) % right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj,obj.ne(1))==0),...
                    jj>obj.ne(1)) % up, left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1-obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj,obj.ne(1))==0),...
                    jj<=(obj.Ne-obj.ne(1))) % up, right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1+obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj-1,obj.ne(1))==0),...
                    jj>obj.ne(1)) % down, left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1-obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj-1,obj.ne(1))==0),...
                    jj<=(obj.Ne-obj.ne(1))) % down, right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1+obj.ne(1);
                vec(ll) = 1;
            end
        end
        
        adj = sparse(ind1,ind2,vec,...
            prod(obj.ne),prod(obj.ne));
        
        if obj.ispartial==1
            adj(obj.missing,:) = [];
            adj(:,obj.missing) = [];
        end
    end
    %=================================================================%
    
    
    
    %== GLOBAL_IDX ===================================================%
    %   Convert 2D grid coordinate to a global coordinate in the grid.
    %   For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %   idx2 is the mobility index.
    function k = global_idx(obj,idx1,idx2)
        k = idx1+(idx2-1)*obj.ne(1);
        
        if obj.ispartial
            b0 = ismember((1:prod(obj.ne))',obj.missing);
            t0 = cumsum(b0); % number of entries missing prior to current index
            
            [~,i0] = intersect(k,obj.missing); % find indices in k that are missing
            
            k = k-t0(k); % reduce indices based on missing indices
            k(i0) = NaN; % return NaN for indices that are missing
        end
    end
    %=================================================================%
    
    
    
    %== TWO_IDX ======================================================%
    %   Convert global grid coordinate to 2D index on grid.
    %   For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %   idx2 is the mobility index.
    function [idx1,idx2] = two_idx(obj,k)
        idx1 = mod(k,obj.ne(1));
        idx1(idx1==0) = obj.ne(1);
        
        idx2 = floor((k-1)./obj.ne(1))+1;
    end
    %=================================================================%
    
    
    
    %== L1 ===========================================================%
    %   Compute the first-order Tikhonov operator.
    %   Form is equiavalent to applying no slope at high-high boundary.
    function [l1] = l1(obj)
        l1 = -diag(sum(tril(obj.adj)))+...
            triu(obj.adj);
        l1(size(obj.adj,1),end) = -1;
            % unity on diagonal in final row for stability in square matrix
            % alternatively, this row can be deleted, however this causes
            % issues during Tikhonov inversion using this method
    end
    %=================================================================%
    
    
    
    %== L2 ===========================================================%
    %   Compute the second-order Tikhonov operator.
    %   Form is equiavalent to applying no slope at grid boundary.
    function [l2] = l2(obj)
        l2 = -diag(sum(obj.adj))+...
            triu(obj.adj)+tril(obj.adj);
    end
    %=================================================================%
    
    
    
    %== GRAD =========================================================%
    %   Calculates the gradient in x. Uses simple first order
    %   differences.
    function out = grad(obj,x)

        [~,dr1,dr2] = obj.dr;

        [grad2,grad1] = gradient(reshape(x,obj.ne));
        grad1 = grad1.*dr1;
        grad2 = grad2.*dr2;

        out = [grad1(:),grad2(:)];

    end
    %=================================================================%
    
    
    
    %== PROJECT ======================================================%
    %   Project x onto current grid. Uses simple linear.
    %   interpolation for this purpose. The parameter 'grid_old'
    %   contains the original grid for the input data x.
    function x = project(obj,grid_old,x)
        
        if grid_old.ispartial==1 % added processing for partial grids
            x = grid_old.partial2full(x);
        end
        
        n1 = length(grid_old.edges{1});
        n2 = length(grid_old.edges{2});
        x = reshape(x,[n1,n2]);
        
        [edges1,edges2] = ndgrid(obj.edges{1},obj.edges{2});
        [edges_old1,edges_old2] = ndgrid(grid_old.edges{1},grid_old.edges{2});
        F = griddedInterpolant(edges_old1,edges_old2,x,'linear','linear');
        
        x = F(edges1,edges2);
        x = x(:);
        
        if obj.ispartial==1 % added processing for partial grids
            x = obj.full2partial(x);
        end
        
    end
    %=================================================================%



    %== TRANSFORM ====================================================%
    %   Function to transform kernel functions. Output is a matrix to 
    %   be multiplied by the original kernel, A.
    function B = transform(obj,grid_old)
        
        for ii=1:obj.dim % loop over both dimensions
            dr_inv = 1./(grid_old.nodes{ii}(2:end)-grid_old.nodes{ii}(1:(end-1)));

            t0 = min(1,max(0,...
                bsxfun(@times,...
                grid_old.nodes{ii}(2:end)-obj.nodes{ii}(1:(end-1))',dr_inv)...
                )); % upper matrix of overlap
            t1 = min(1,max(0,...
                bsxfun(@times,...
                obj.nodes{ii}(2:end)'-grid_old.nodes{ii}(1:(end-1)),dr_inv)...
                )); % lower matrix of overlap
            t2{ii} = sparse(min(t0,t1));
        end % end loop over both dimensions

        ind_old_tot = 1:prod(grid_old.ne);
        ind1_old = mod(ind_old_tot-1,grid_old.ne(1))+1;
        ind2_old = ceil(ind_old_tot./grid_old.ne(1));
        
        ind_tot = 1:prod(obj.ne);
        ind1 = mod(ind_tot-1,obj.ne(1))+1;
        ind2 = ceil(ind_tot./obj.ne(1));
        
        B = (t2{1}(ind1,ind1_old).*t2{2}(ind2,ind2_old))';
        
        %-- Added processing for partial grids -----------------------%
        if obj.ispartial==1; B(:,obj.missing) = []; end
        if grid_old.ispartial==1; B(grid_old.missing,:) = []; end
        %-------------------------------------------------------------%
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
        dr1 = abs(dr1); % in case edges vector is reversed
        dr2 = abs(dr2);
        dr = dr1(:).*dr2(:);
        
        %-- Added processing for partial grids -----------------------%
        if obj.ispartial==1
            [~,r_min,r_max] = obj.ray_sum([0,obj.cut(1)],obj.cut(2),0);
            t0 = (r_min(:,1)-log10(obj.nelements(:,1))).*...
                (log10(obj.nelements(:,2))-log10(obj.nelements(:,1)));
                    % lower rectangle
            t1 = (log10(obj.nelements(:,4))-r_max(:,2)).*...
                (log10(obj.nelements(:,2))-r_min(:,1));
                    % right rectangle
            t2 = 1/2.*(r_max(:,1)-r_min(:,1)).*...
                (r_max(:,2)-r_min(:,2));
                    % upper, left triangle
            dr = t0+t1+t2;
            
            if length(obj.cut)==4 % Note: Ignores if both rays pass through element
                [~,r_min,r_max] = obj.ray_sum([0,obj.cut(3)],obj.cut(4),0);
                t0 = (log10(obj.nelements(:,2))-r_max(:,1)).*...
                    (log10(obj.nelements(:,2))-log10(obj.nelements(:,1)));
                        % upper rectangle
                t1 = (r_min(:,2)-log10(obj.nelements(:,3))).*...
                    (r_max(:,1)-log10(obj.nelements(:,1)));
                        % left rectangle
                t2 = 1/2.*(r_max(:,1)-r_min(:,1)).*...
                    (r_max(:,2)-r_min(:,2));
                        % lower, right triangle
                dr = t0+t1+t2;
            end
            
        end
        %-------------------------------------------------------------%
    end
    %=================================================================%


    
    %== MARGINALIZE ==================================================%
    %   Marginalizes over the grid in each dimension.
    %   Uses Euler's method to integrate over domain.
    function [marg,tot] = marginalize(obj,x)
        
        x = obj.reshape(x);
        
        [dr,dr1,dr2] = obj.dr; % generate differential area of elements
        dr = obj.reshape(dr); % fills out dr for partial grids
        
        tot = sum(x(:).*dr(:)); % integrated total
        
        marg{1} = sum(dr2.*x,2); % integrate over diameter
        marg{2} = sum(dr1.*x,1); % integrate over mass
        
    end
    %=================================================================%
    
    
    
    %== RESHAPE ======================================================%
    %   A simple function to reshape a vector based on the grid.
    %   Note: If the grid is partial, missing grid points are 
    %   filled with zeros. 
    function x = reshape(obj,x)
        
        if obj.ispartial==1 % if partial grid
            x = obj.partial2full(x);
        end
        
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
        
        if nargout>1; t1 = obj.elements(:,1); end
        if nargout>2; t2 = obj.elements(:,2); end
    end
    %=================================================================%



    %== RAY_SUM ======================================================%
    %   Perfrom a ray sum for a given ray and the current grid.
    %   Currently assumes uniform, logarithmic grid 
    %   and can accommodate partial grids.
    %   Based on:	Code from Samuel Grauer
    %   Author:     Timothy Sipkens, 2019-07-14
    %-----------------------------------------------------------------%
    % Inputs:
    %   logr0   A single point on the line in log-log space, r0 = log10([dim1,dim2])
    %	slope   Slope of the line
    %   f_bar   Flag for progress bar
    % Outputs:
    %   C       Ray-sum matrix
    %-----------------------------------------------------------------%
    function [C,rmin,rmax] = ray_sum(obj,logr0,slope,f_bar)
        
        if ~exist('f_bar','var'); f_bar = []; end
        if isempty(f_bar); f_bar = 0; end
        
        %-- Preallocate arrays ---------------------------------------%
        m = size(slope,1);
        C = spalloc(m,obj.Ne,ceil(0.1*m*obj.Ne)); % assume 10% full
        
        
        %-- Compute ray-sum matrix -----------------------------------%
        if f_bar; tools.textbar(); end
        for ii=1:m % loop over multiple rays

            %-- Ray vector -------------%
            dv = [slope,1]; % convert slope to step vector along line
            dv = dv/norm(dv);
            dv(dv == 0) = 1e-10; % for stability during division
            
            
            %-- Line intersections -----%
            %   Use parametric representation of the line and find two
            %   intersections or each element.
            %   Note: Assumes a logarithmic grid.
            ttmp = (log10(obj.nelements(:,[1,3]))-...
                logr0)./dv; % minimum of element
            tmax = (log10(obj.nelements(:,[2,4]))-...
                logr0)./dv; % maximum of element
            
            
            %-- Corrections ------------%
            %   Decide which points would correspond to transecting the
            %   pixel.
            tmin = max(min(ttmp,tmax),[],2);
            tmax = min(max(ttmp,tmax),[],2);
            tmax(tmax<tmin) = tmin(tmax<tmin); % check if crosses pixel
            
            
            %-- Convert back to [x,y] --%
            rmin  = logr0+tmin.*[dv(2),dv(1)]; % location of intersect with min. of pixel
            rmax = logr0+tmax.*[dv(2),dv(1)]; % location of intersect with max. of pixel
            chord = sqrt(sum((rmax-rmin).^2,2)); % chord length
            chord(chord<1e-15) = 0; % truncate small values
            
            
            %-- Ray-sum matrix ---------%
            [~,jj,a] = find(chord');
            if ~isempty(a)
                C(ii,:) = sparse(1,jj,a,1,obj.Ne,ceil(0.6*obj.Ne));
            end
            if f_bar, tools.textbar(ii/m); end

        end % end loop over multiple rays

    end
    %=================================================================%
    
    
    
    %== CLOSEST_IDX =================================================%
    %   Returns the pixel in which r0 is located.
    %   This function uses vector operations to find multiple points.
    %-----------------------------------------------------------------%
    % Inputs:
    %   r0      Coordinates in grid space, r0 = [dim1,dim2]
    %           Can form N x 2 vector, where N is the number of points
    %           to be found.
    % Outputs:
    %   idx     Global index on the grid, incorporating missing pixels
    %   idx_2d  Pair of indices of pixel location
    %-----------------------------------------------------------------%
    function [k,idx_2d] = closest_idx(obj,r0)
        
        idx_2d = zeros(size(r0,1),2); % pre-allocate
        
        for ii=1:obj.dim
            idx_bool = and(...
                r0(:,ii)>=obj.nodes{ii}(1:(end-1)),...
                r0(:,ii)<obj.nodes{ii}(2:end));
                % if on border, place in higher pixel
            
            [~,idx_2d(:,ii)] = max(idx_bool,[],2);
            
            idx_2d(~any(idx_bool,2),ii) = NaN;
                % if point is outside of grid, return NaN
        end
        
        f_nan = or(isnan(idx_2d(:,1)),isnan(idx_2d(:,2)));
        k = NaN(size(r0,1),1);
        k(~f_nan) = obj.global_idx(idx_2d(~f_nan,1),idx_2d(~f_nan,2));
    end
    %=================================================================%

    
    
%=====================================================================%
%-- VISUALIZATION METHODS --------------------------------------------%
%=====================================================================%
    
    %== PLOT2D =======================================================%
    %   Plots x as a 2D function on the grid.
    %   Author: Timothy Sipkens, 2018-11-21
    function [h,x] = plot2d(obj,x,f_contf)
        
        if ~exist('f_contf','var'); f_contf = []; end % set empty contourf flag
        if isempty(f_contf); f_contf = 0; end % set contourf flag to false
        
        %-- Issue warning if grid edges are not be uniform -----------%
        %   The imagesc function used here does not conserve proportions.
        [~,dr1,dr2] = obj.dr; % used to give warning below
        dr0 = dr1(:).*dr2(:);
        if ~all(abs(dr0(2:end)-dr0(1))<1e-10)
            warning(['The plot2d method does not display ',...
                'correct proportions for non-uniform grids.']);
        end
        
        x = obj.reshape(x);
        
        %-- Plot -------------------------------%
        if f_contf==0 % plot as image
            imagesc(obj.edges{2},obj.edges{1},x);
            set(gca,'YDir','normal');
        else % plot as contourf
            contourf(obj.edges{2},obj.edges{1},x,15,'EdgeColor','none');
        end
        
        %-- Adjust tick marks for log scale ----%
        if strcmp('logarithmic',obj.discrete)
            set(gca,'XScale','log');
            set(gca,'YScale','log');
        end
        
        xlim(obj.span(2,:));
        ylim(obj.span(1,:));
        
        if nargout>0; h = gca; end
        
    end
    %=================================================================%



    %== PLOT2D_MARG ==================================================%
    %   Plots x as a 2D function on the grid, with marginalized distributions.
    %   Author: Timothy Sipkens, 2018-11-21
    function [h,x_m] = plot2d_marg(obj,x,obj_t,x_t,f_contf)
        
        if ~exist('f_contf','var'); f_contf = []; end % set empty contourf flag
        if ~exist('x_t','var'); x_t = []; end
        
        subplot(4,4,[5,15]);
        obj.plot2d(x,f_contf);
        
        x_m = obj.marginalize(x);
        
        
        %-- Plot marginal distribution (dim 2) -----------------------%
        subplot(4,4,[1,3]);
        marg_dim = 2;
        stairs(obj.nodes{marg_dim},...
            [x_m{marg_dim},0],'k');
        xlim([min(obj.edges{marg_dim}),max(obj.edges{marg_dim})]);
        set(gca,'XScale','log');
        
        %-- Also plot marginal of the true distribution --------------%
        if ~isempty(x_t)
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
        
        %-- Also plot marginal of the true distribution --------------%
        if ~isempty(x_t)
            hold on;
            plot([0;x_m_t{marg_dim}],...
                obj_t.nodes{marg_dim},'color',[0.6,0.6,0.6]);
            hold off;
        end
        
        subplot(4,4,[5,15]);
        if nargout>0; h = gca; end
        
    end
    %=================================================================%



    %== PLOT2D_SWEEP =================================================%
    %   Plot data in slices, sweeping through the provided colormap.
    %   Author: Timothy Sipkens, 2019-11-28
    function [h,x] = plot2d_sweep(grid,x,cmap)
        
        n1 = ceil(grid.ne(1)./20);
        n2 = floor(grid.ne(1)/n1);
        n3 = floor(length(cmap)/n2);
        cmap2 = cmap(1:n3:end,:);

        set(gca,'ColorOrder',cmap2,'NextPlot','replacechildren');
        x = reshape(x,grid.ne);
        h = semilogx(grid.edges{2},x(1:n1:end,:),...
            'o-','MarkerSize',2.5,'MarkerFaceColor',[1,1,1]);

        if nargout==0; clear h; end
    
    end
    %=================================================================%
    
    
    
    %== PLOT2D_SWEEPT ================================================%
    %   Plot transposed data in slices, sweeping through the provided colormap.
    %   Author: Timothy Sipkens, 2019-11-28
    function [h,x] = plot2d_sweept(grid,x,cmap)
        
        n1 = ceil(grid.ne(2)./20);
        n2 = floor(grid.ne(2)/n1);
        n3 = floor(length(cmap)/n2);
        cmap2 = cmap(1:n3:end,:);
        
        set(gca,'ColorOrder',cmap2,'NextPlot','replacechildren');
        x = reshape(x,grid.ne)';
        h = semilogx(grid.edges{1},x(1:n1:end,:),...
            'o-','MarkerSize',2.5,'MarkerFaceColor',[1,1,1]);
        
        if nargout==0; clear h; end
        
    end
    %=================================================================%
    
    
    
    %== PLOT_MARGINAL ================================================%
    %   Plot marginal distributions
    %   Author:	Timothy Sipkens, 2019-07-17
    %   Note: 'x' can be a cell array containing multiple x vectors
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
    
    
    
    
%=====================================================================%
%-- SUPPORT FOR PARTIAL GRIDS ----------------------------------------%
%=====================================================================%
    
    %== PARTIAL ======================================================%
    %   Convert to a partial grid. Currently takes a y-intercept, r0, 
    %   and slope0 as arguements and cuts upper triangle.
    %   Added r1 and slope1 arguments will also cut a lower triangle.
    function obj = partial(obj,r0,slope0,r1,slope1)
        
        %-- Parse inputs ----------------------------%
        if ~exist('slope0','var'); slope0 = []; end
        if isempty(slope0); slope0 = 1; end
        
        if ~exist('r0','var'); r0 = []; end
        if length(r0)==1; b0 = r0; end % if scalar, use as y-intercept
        if length(r0)==2; b0 = r0(1)-slope0*r0(2); end % if coordinates, find y-intercept
        if isempty(r0); b0 = 0; end % if not specified, use b = 0
        
        %-- For bottom triangle --%
        if ~exist('slope1','var'); slope1 = []; end
        if isempty(slope1); slope1 = 0; end
        
        if ~exist('r1','var'); r1 = -inf; end
        if length(r1)==1; b1 = r1; end % if scalar, use as y-intercept
        if length(r1)==2; b1 = r1(1)-slope1*r1(2); end % if coordinates, find y-intercept
        if isempty(r1); b1 = 0; end % if not specified, use b = 0
        %--------------------------------------------%
        
        
        if strcmp(obj.discrete,'logarithmic')
            t0 = log10(obj.elements);
        else
            t0 = obj.elements;
        end
        
        %-- Cut upper triangle ---------------------%
        f_missing = t0(:,1)>(t0(:,2).*slope0+b0);
        t1 = 1:prod(obj.ne);
        obj.cut = [b0,slope0];
        
        
        %-- Consider cutting lower triangle --------%
        if ~isinf(r1)
            f_missing1 = t0(:,1)<(t0(:,2).*slope1+b1);
            f_missing = or(f_missing,f_missing1);
            obj.cut = [obj.cut,b1,slope1];
        end
        
        
        %-- Update grid properties -----------------%
        obj.ispartial = 1;
        obj.missing = t1(f_missing);
        
        obj.elements = obj.elements(~f_missing,:);
        obj.nelements = obj.nelements(~f_missing,:);
        obj.Ne = size(obj.elements,1);
        obj = obj.padjacency;
        
    end
    %=================================================================%
    
    
    
    %== PARTIAL2FULL =================================================%
    %   Convert x defined on a partial grid to the full grid equivalent, 
    %   using zeros to fill the removed grid points.
    function x_full = partial2full(obj,x)
        x_full = zeros(prod(obj.ne),1);
        t0 = setdiff((1:prod(obj.ne))',obj.missing);
        x_full(t0) = x;
    end
    %=================================================================%
    
    
    
    %== FULL2PARTIAL =================================================%
    %   Convert x defined on a full grid to the partial grid equivalent, 
    %   removing entries for missing indices.
    function x = full2partial(obj,x)
        x(obj.missing,:) = [];
    end
    %=================================================================%
    
    
    
    %== PADJACENCY ===================================================%
    %   Convert x defined on a partial grid to the full grid equivalent, 
    %   using zeros to fill the removed grid points.
    function [obj,adj] = padjacency(obj)
        [~,adj] = obj.adjacency;
        
        adj(obj.missing,:) = [];
        adj(:,obj.missing) = [];
        
        obj.adj = adj;
    end
    %=================================================================%
    
end

end
