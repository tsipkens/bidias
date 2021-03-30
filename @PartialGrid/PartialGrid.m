
% PARTIALGRID  A subclass of Grid handling partially truncated grids.
%  
%  G = PartialGrid(...,R0,SLOPE0) tkaes the same inputs as the Grid
%  class constructor, but adds extra inputs to form a partial grid. 
%  This call removes elements above the line that goes through the 
%  point R0 and having a slope of SLOPE0. For logarithmic
%  grids, lines correspond to exponential curves, R0 are given as 
%  log10(...) quantities, and slopes correspond to the exponent. For 
%  example, Grid.partial([0,2],3) removes all grid elements above the
%  exponential curve that passes through [1,100] and having an exponent
%  of 3 (e.g., curves that increase volumetrically). 
%  
%  G = PartialGrid(...,R0,SLOPE0,R1,SLOPE1) adds a second set of 
%  arguments analogous to above but remove points below a given line 
%  (instead of above). 
%  
%  AUTHOR: Timothy Sipkens, 2021-03-29

classdef PartialGrid < Grid
   
    
    
properties
    
    missing = [];  % global indices of missing grid points for partial grids
    cut = [];      % [y-intercept,slope] used to cut partial grid
    
end
    


methods
        
    %== PARTIALGRID ======================================================%
    %   Constructor class.
    function obj = PartialGrid(span_edges, ne, discrete, r0, slope0, r1, slope1)
        
        % Call Grid constructor.
        obj = obj@Grid(span_edges, ne, discrete);

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

        %-- Cut upper triangle ---------------------%
        if strcmp(obj.discrete,'log')
            tup = log10(obj.nelements(:,[1,4]));
        else
            tup = obj.nelements(:,[1,4]);
        end
        tup = tup+abs(1e-3.*mean(tup(2:end,:)-tup(1:(end-1),:))).*[1,0];
            % avoids minimially overlapping elements

        f_missing = tup(:,1) > (tup(:,2).*slope0 + b0);
        t1 = 1:prod(obj.ne);
        obj.cut = [b0, slope0];


        %-- Consider cutting lower triangle --------%
        if strcmp(obj.discrete,'log')
            tlow = log10(obj.nelements(:,[2,3]));
        else
            tlow = obj.nelements(:,[2,3]);
        end
        tlow = tlow+abs(1e-3.*mean(tlow(2:end,:)-tlow(1:(end-1),:))).*[-1,0];
            % avoids minimially overlapping elements

        if ~isinf(r1)
            f_missing1 = tlow(:,1)<(tlow(:,2).*slope1+b1);
            f_missing = or(f_missing,f_missing1);
            obj.cut = [obj.cut,b1,slope1];
        end


        %-- Update grid properties -----------------%
        obj.missing = t1(f_missing);

        obj.elements = obj.elements(~f_missing,:);
        obj.nelements = obj.nelements(~f_missing,:);
        obj.Ne = size(obj.elements,1);
        obj = obj.adjacency;

    end
    
    
    
    %== ADJACENCY ====================================================%
    %   Compute the adjacency matrix for the full grid,
    %   using a four-point stencil.
    function [obj,adj] = adjacency(obj)
        
        [~, adj] = adjacency@Grid(obj);
        
        adj(obj.missing, :) = [];
        adj(:, obj.missing) = [];
        
        obj.adj = adj;
        
    end
    
    
    
    %== GLOBAL_IDX ===================================================%
    %   Convert 2D grid coordinate to a global coordinate in the grid.
    %   For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %   idx2 is the mobility index.
    function k = global_idx(obj, idx1, idx2)
        
        k = global_idx@Grid(obj, idx1, idx2);  % use Grid function
        
        %-- Update for partial grid. --%
        b0 = ismember((1:prod(obj.ne))',obj.missing);
        t0 = cumsum(b0); % number of entries missing prior to current index

        [~,i0] = intersect(k,obj.missing); % find indices in k that are missing

        k = k-t0(k); % reduce indices based on missing indices
        k(i0) = NaN; % return NaN for indices that are missing
    end
    %=================================================================%
    
    
    
    %== PROJECT ======================================================%
    %   Project x onto current grid. Uses simple linear.
    %   interpolation for this purpose. The parameter 'grid_old'
    %   contains the original grid for the input data x.
    function x = project(obj, grid_old, x)
        
        x = project@Grid(obj, grid_old, x);
        
        x = obj.full2partial(x);
        
    end
    %=================================================================%
    
    
    
    %== TRANSFORM ====================================================%
    %   Function to transform kernel functions. Output is a matrix to 
    %   be multiplied by the original kernel, A.
    function B = transform(obj, grid_old)
        
        B = transform@Grid(obj, grid_old);
        
        B(:,obj.missing) = [];
        if isa(grid_old, 'PartialGrid')
            B(grid_old.missing,:) = [];
        end
        %-------------------------------------------------------------%
    end
    %=================================================================%
    
    
    
    %== DR ===========================================================%
    %   Calculates the differential area of the elements in the grid.
    function [dr,dr1,dr2] = dr(obj)
        
        [dr,dr1,dr2] = dr@Grid(obj);
        
        %-- Added processing for partial grids -----------------------%
        dr0 = obj.full2partial(dr); % used if lower cut is employed

        [~,rmin,rmax] = obj.ray_sum([0,obj.cut(1)],obj.cut(2),0);
        t0 = (rmin(:,1)-log10(obj.nelements(:,1))).*...
            (log10(obj.nelements(:,4))-log10(obj.nelements(:,3)));
                % lower rectangle
        t1 = (log10(obj.nelements(:,4))-rmax(:,2)).*...
            (log10(obj.nelements(:,2))-rmin(:,1));
                % right rectangle
        t2 = 1/2.*(rmax(:,1)-rmin(:,1)).*...
            (rmax(:,2)-rmin(:,2));
                % upper, left triangle
        dr = t0+t1+t2;

        if length(obj.cut)==4 % consider cutting lower triangle
            [~,rmin,rmax] = obj.ray_sum([0,obj.cut(3)],obj.cut(4),0);
            t0 = (log10(obj.nelements(:,2))-rmax(:,1)).*...
                (log10(obj.nelements(:,4))-log10(obj.nelements(:,3)));
                    % upper rectangle
            t1 = (rmin(:,2)-log10(obj.nelements(:,3))).*...
                (rmax(:,1)-log10(obj.nelements(:,1)));
                    % left rectangle
            t2 = 1/2.*(rmax(:,1)-rmin(:,1)).*...
                (rmax(:,2)-rmin(:,2));
                    % lower, right triangle
            dr = dr.*(t0+t1+t2)./dr0; % accounts for element that are discected twice
        end
        %-------------------------------------------------------------%
    end
    %=================================================================%
    
    
    
    %== MARGINALIZE_OP ===============================================%
    %   A marginalizing operator, C1, to act on 2D distributions.
    %   Author: Timothy Sipkens, Arash Naseri, 2020-03-09
    function [C1,dr0] = marginalize_op(obj, dim)
        
        C1 = marginalize_op@Grid(obj, dim);
        
        C1(:,obj.missing) = []; % remove missing elements

        %-- Extra processing to account for partial elements -----%
        [dr,dr1,dr2] = obj.dr;
        switch dim % determine which dimension to sum over
            case 2 % integrate in column direction
                dr1 = obj.full2partial(dr1(:));
                dr0 = dr./dr1; % ~dr2
            case 1 % integrate in row direction
                dr2 = obj.full2partial(dr2(:));
                dr0 = dr./dr2; % ~dr1
        end
        C1 = bsxfun(@times,C1,dr0');
        
    end
    %=====================================================================%
    
    
    
    %== RESHAPE ======================================================%
    %   A simple function to reshape a vector based on the grid.
    %   Note: If the grid is partial, missing grid points are 
    %   filled with zeros. 
    function x = reshape(obj, x)
        
        x = obj.partial2full(x);
        
        x = reshape@Grid(obj, x);
        
    end
    %=================================================================%
    
    
    
    %== PLOT2D =======================================================%
    %   Plots x as a 2D function on the grid.
    %   Author: Timothy Sipkens, 2018-11-21
    function [h, x] = plot2d(obj, x, f_contf)
        
        if ~exist('f_contf', 'var'); f_contf = []; end
        
        % Use default method to start plot.
        [h, x] = plot2d@Grid(obj, x, f_contf);
        
        % Add lines marking the edges of the partial grid.
        hold on;
        tools.overlay_line(obj, [0,obj.cut(1)], obj.cut(2), ...
            'Color', [0.5,0.5,0.5]); % overlay partial grid limits, gray lines

         % If also a bottom cut.
        if length(obj.cut)>2
            tools.overlay_line(obj, [0,obj.cut(3)], obj.cut(4), ...
                'Color', [0.5,0.5,0.5]); % add a gray line
        end
        hold off;
        
        % Reapply grey labels and axes to allow viz against dark and light bgs.
        set(gca, 'XColor', [0.5, 0.5, 0.5], ...
            'YColor', [0.5, 0.5, 0.5], ...
            'linewidth', 0.75);
        
        if nargout>0; h = gca; end
        
    end
    %=================================================================%

    
    
    %== PARTIAL2FULL =================================================%
    %   Convert x defined on a partial grid to the full grid equivalent, 
    %   using zeros to fill the removed grid points.
    function x_full = partial2full(obj ,x)
        x_full = zeros(prod(obj.ne),1);
        t0 = setdiff((1:prod(obj.ne))', obj.missing);
        x_full(t0) = x;
    end
    %=================================================================%
    
    
    
    %== FULL2PARTIAL =================================================%
    %   Convert x defined on a full grid to the partial grid equivalent, 
    %   removing entries for missing indices.
    function x = full2partial(obj, x)
        x(obj.missing,:) = [];
    end
    %=================================================================%
    
end

end

