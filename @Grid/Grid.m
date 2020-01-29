
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
    Ne = [];    % total number of pixels/elements, i.e. prod(ne)
    
    edges = []; % cell of vectors containing edge points of pixel/element
                % centers, with one cell entry per dimension
    
    elements = []; % contains position of pixel/element centers as a (ne x 2) vector
    nodes = []; % contains position of nodes surrounding elements for each dimension
                % each cell has a vector of size (ne + 1).
                
    adj = [];   % adjacency matrix
    
    ispartial = 0; % toggle of whether the grid is 'partial' or sparse,
                   % that is having some grid points missing
    missing = [];  % global indices of missing grid points for partial grids
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
        
        if nargin==0; return; end % return empty grid

        if isa(span_edges,'cell') % consider case where edges are given
            obj.edges = span_edges;
            obj.ne = [length(span_edges{1}),...
                length(span_edges{2})];
            obj.span = [min(span_edges{1}),max(span_edges{1});...
                min(span_edges{2}),max(span_edges{2})];
            obj.discrete = 'custom';

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

        %-- Generate vectorized list of elements ---------------------%
        %   One column per dimension
        [vec1{1},vec1{2}] = ndgrid(obj.edges{1},obj.edges{2});
        obj.elements(:,1) = vec1{1}(:); % vectorize output
        obj.elements(:,2) = vec1{2}(:);
    end
    %=================================================================%
    
    
    
    %== ADJACENCY ====================================================%
    %   Compute the adjacency matrix for the full grid.
    function [obj,adj] = adjacency(obj)
        
        ind1 = zeros(3*prod(obj.ne),1);
        ind2 = ones(3*prod(obj.ne),1);
        vec = zeros(3*prod(obj.ne),1);
        ll = 0;
        
        for jj=1:prod(obj.ne)
            if ~(mod(jj,obj.ne(1))==0)
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1;
                vec(ll) = 1;
            end
            
            if ~(mod(jj-1,obj.ne(1))==0)
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1;
                vec(ll) = 1;
            end
            
            if jj>obj.ne(1)
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-obj.ne(1);
                vec(ll) = 1;
            end
            
            if jj<=(obj.Ne-obj.ne(1))
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
    
    
    
    %== GLOBAL_IDX ===================================================%
    %   Convert 2D grid coordinate to a global coordinate in the grid.
    %   For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %   idx2 is the mobility index.
    function k = global_idx(obj,idx1,idx2)
        k = idx1+(idx2-1)*obj.ne(1);
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
    %   Function to integrate over uniform basis functions to transform 
    %   kernel functions. Output is a matrix to be multiplied by original
    %   kernel, A.
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
        dr = dr2(:).*dr1(:);
        
        dr = obj.full2partial(dr); % added processing for partial grids
        
    end
    %=================================================================%


    
    %== MARGINALIZE ==================================================%
    %   Marginalizes over the grid in each dimension.
    %   Uses Euler's method to integrate over domain.
    function [marg,tot] = marginalize(obj,x)
        
        x = obj.reshape(x);
        
        [~,dr1,dr2] = obj.dr; % generate differential area of elements
        dr = dr2(:).*dr1(:);
        
        tot = sum(x(:).*dr); % integrated total
        
        marg{1} = sum(dr2.*x,2); % integrate over diameter
        marg{2} = sum(dr1.*x,1); % integrate over mass
        
    end
    %=================================================================%
    
    
    
    %== RESHAPE ======================================================%
    %   A simple function to reshape a vector based on the grid.
    %   Note: If the grid is partial, missing grid points are 
    %           filled with zeros. 
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
            tmax(tmax<tmin) = tmin(tmax<tmin); % check if crosses pixel


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

    
    
%=====================================================================%
%-- VISUALIZATION METHODS --------------------------------------------%
%=====================================================================%
    
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
    
    
    
%=====================================================================%
%-- SUPPORT FOR PARTIAL GRIDS ----------------------------------------%
%=====================================================================%
    
    %== PARTIAL =====================================================%
    %   Convert to a partial grid. Currently take a y-intercept, b, 
    %   and slope, m, as arguements and cuts upper triangle.
    function obj = partial(obj,b,m)
        
        if ~exist('b','var'); b = []; end
        if isempty(b); b = 0; end
        
        if ~exist('m','var'); m = []; end
        if isempty(m); m = 0; end
        
        if strcmp(obj.discrete,'logarithmic')
            t0 = log10(obj.elements);
        else
            t0 = obj.elements;
        end
        
        bool_above = t0(:,1)>(t0(:,2).*m+b);
        t1 = 1:length(bool_above);
        
        %-- Update grid properties -----------------%
        obj.ispartial = 1;
        obj.missing = t1(bool_above);
        obj.elements = obj.elements(~bool_above,:);
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
        idx_miss = sort(obj.missing,'descend');
        x(idx_miss,:) = [];
    end
    %=================================================================%
    
    
    
    %== PADJACENCY ===================================================%
    %   Convert x defined on a partial grid to the full grid equivalent, 
    %   using zeros to fill the removed grid points.
    function [obj,adj] = padjacency(obj)
        [~,adj] = obj.adjacency;
        
        idx_miss = sort(obj.missing,'descend');
        
        adj(idx_miss,:) = [];
        adj(:,idx_miss) = [];
        
        obj.adj = adj;
    end
    %=================================================================%
    
end


methods(Static)
    %== CLOSEST_IDX =================================================%
    %   Returns the effective index closest to r position.
    function idx = closest_idx(edges,r)

        dim = length(edges); % will be number of dimensions

        for ii=1:dim
            idx_bool = and(...
                r(:,ii)>=edges{ii}(1:(end-1)),...
                r(:,ii)<edges{ii}(2:end));
            
            idx_temp = sum(cumprod(idx_bool==0,2),2)+1;
            idx_temp(r(:,ii)==edges{ii}(end)) = ...
                idx_temp(r(:,ii)==edges{ii}(end))-1;
                    % condition for points at end of domain
            idx(:,ii) = idx_temp;
        end
    end
    %=================================================================%
end
end
