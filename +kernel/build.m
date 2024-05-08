
% BUILD  A general function to compute a 2D transfer function. 
%  
%  NOTE: This function generally requires that the type of particle size is
%  specified in grid.type = {}. The default if not specified is to assign
%  grid.type = {'dm', 'mp'}, corresponding to a mass-mobility scenario.
%  Other cases require identifying the dimensions of the grid. 
%  
%  A = kernel.build_grid(grid_i, z_vec, ...) builds a kernel based
%  on the grid for evaluation (GRID_I). Setpoints should be specified 
%  for each classifer Charge states to evaluate at are provided in Z_VEC. 
%  When Z_VEC is empty, the function defaults to Z_VEC = 1:3. Remaining 
%  elements are required and consistitute name-cell pairs for each classifier. 
%  
%  [A,AC] = kernel.build_grid(...) adds an output AC that is not summed
%  ovre the charge states, which is useful for advanced analysis. 
%  
%  A = kernel.build_grid(..., 'pma', {M_STAR, PROP_P, VARARGIN})
%  Builds a PMA contribution to the kernel using the properties in PROP_P
%  and a set of secondary information for the setpoints 
%  (e.g., {M_STAR, PROP_P, 'Rm', Rm}). M_STAR is expected in fg. 
%  
%  A = kernel.build_grid(..., 'dma', {D_STAR, PROP_D, VARARGIN})
%  Builds a DMA contribution to the kernel using the properties in PROP_D
%  and a set of secondary information for the classifer evaluations in
%  VARARGIN that is passed directly to tfer_dma. D_STAR is expected in nm.
%  
%  A = kernel.build_grid(..., 'charger', {VARARGIN})
%  Builds a charger contribution to the transfer function. An empty cell
%  will use the default call to charger(...). This addition is required for
%  all classifiers that require charging for classification. 
%  This does not require a mobility diameter setpoint. 
%  
%  Other classifiers follow this template. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2024-01-11

function [A, Ac] = build(grid_i, z_vec, varargin)

if mod(length(varargin), 2) ~= 0; error('Wrong number of inputs.'); end

if ~exist('z_vec', 'var'); z_vec = []; end
if isempty(z_vec); z_vec = 1:3; end  % default charge states for evaluation

tools.textheader('Computing kernel')

nc = (length(varargin)/2);  % number of classifiers
Lambda{nc} = [];  % initialize Lambda_ii

dm_idx = find(strcmp(grid_i.type, 'dm'));
mp_idx = find(strcmp(grid_i.type, 'mp'));
da_idx = find(strcmp(grid_i.type, 'da'));

% Set default indices, if type not specified (mass-mobility grid).
if isempty(grid_i.type)
    dm_idx = 2;
    mp_idx = 1;
end

% Handle if mobility diameter is not given directly (req'd for charging/PMA).
% Compute using known relationships.
if isempty(dm_idx)
    addpath autils;  % add aerosol utilities (autils) package
    
    % OPTION 1: Use da and mp to directly compute dm.
    if and(~isempty(da_idx), ~isempty(mp_idx))
        m = grid_i.elements(:, mp_idx);  % get mass from relevant dimension
        da = grid_i.elements(:, da_idx);  % get da from relevant dimension
        dm = mp_da2dm(m, da);  % fully constrained calculation
        dm2 = dm';
    
    % OPTION 2: Apply assumption of a mass-mobility relationship, 
    % which will be less precise. Necessary for PMA-SP2.
    elseif ~isempty(mp_idx)
        idx_p = find(strcmp(varargin, 'pma')) + 1;  % first find PMA inputs
        prop_p = varargin{idx_p}{2};  % extract prop_pma from PMA input
        
        % Then convert using mass-mobility relationship. 
        disp(' Invoking mass-mobility relationship to determine dm.');
        m = grid_i.elements(:, mp_idx);  % get mass from relevant dimension
        dm = mp2dm(m .* 1e-18, prop_p) .* 1e9;
        dm2 = dm';
    end
end


% Loop over the various classifiers. 
for ii=1:nc
    jj = 2*ii - 1;  % index of classifier

    switch varargin{jj}
        
        %== CHARGER ======================================================%
        case 'charger'
            disp(' Computing charger contribution ...');
            
            % Assign inputs.
            if isempty(dm_idx)  % then PMA without DMA
                dm = dm';  % inherit from previous PMA calc.
                dm2 = dm';
            else
                dm = grid_i.edges{dm_idx};
                dm2 = grid_i.elements(:, dm_idx);
            end

            f_z = charger(dm', z_vec, varargin{jj+1}{1:end}); % get fraction charged for d vector
            Lambda{ii} = permute(f_z, [3, 2, 1]);

            % Duplicate over other grid dimensions.
            [~,kk] = max(dm == dm2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            tools.textdone();


        %== SMPS =========================================================%
        case {'dma', 'smps'}
            disp(' Computing DMA contribution ...');
            
            % Unpack inputs.
            d_star2 = varargin{jj+1}{1};
            d_star = unique(d_star2)';
            d = grid_i.edges{dm_idx};  % points for integration
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_dma(d_star, d', z_vec, varargin{jj+1}{2:end});

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, dm_idx);
            [~,kk] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            [~,kk] = max(d_star == d_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();

        
        %== PMA ==========================================================%
        %   Currently assumes other dimension is mobility diameter. 
        case 'pma'
            disp(' Computing PMA contribution ...');

            % Unpack inputs.
            m_star = varargin{jj+1}{1};  % don't use unique(), as resolution may change
            prop_p = varargin{jj+1}{2};  % DMA properties

            % Points for integration.
            m = grid_i.elements(:, mp_idx);  % masses at which to compute the transfer function (not setpoints)

            % Handle mobility diameter.
            if ~isempty(dm_idx)  % use corresponding dimension of grid
                dm = grid_i.elements(:, dm_idx);
            elseif length(varargin{jj+1}) > 2
                if ~isempty(varargin{jj+1}{3})
                    dm = varargin{jj+1}{3};  % then get from explicit input
                end
            else  % then likely PMA without DMA
                dm = (m .* 1e-18 ./ prop_p.rho0) .^ ...
                    (1/prop_p.Dm) .* 1e9;  % use mass-mobility
                disp('  Invoking mass-mobility relation for PMA.')
            end

            addpath 'tfer\tfer-pma';  % added to calculate sp
            sp = get_setpoint(prop_p,...  % get PMA setpoints
                'm_star', m_star .* 1e-18, ...  % mass from the grid
                varargin{jj+1}{3:end});  % extra name-value pair to specify setpoint
            
            % Find unique setpoints. 
            [spu, ~, kk] = unique([[sp.m_star]', [sp.Rm]'], 'rows');
            spu = get_setpoint(prop_p, 'm_star', spu(:,1), 'Rm', spu(:,2));
            
            Lambda{ii} = tfer_pma(...
                spu, m, dm, z_vec, prop_p);

            % Duplicate over repeat entries.
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();
            

        %== AAC ==========================================================%
        case 'aac'
            disp(' Computing AAC contribution ...');
            
            addpath 'autils';
            
            % Assign inputs.
            d_star2 = varargin{jj+1}{1};

            % Copy over prop structure.
            prop = varargin{jj+1}{2};
            f = {'Qa', 'Qs', 'Qsh', 'Qexh'};
            nf = length(f);
            
            % Get unique flow/da_star combinations.
            % Build unique vector.
            u = d_star2;
            for ff=1:nf  % extend flows to allow for vector input
                u = [u, prop.(f{ff}) .* ones(size(d_star2))];
            end
            [u, ~, kk] = unique(u, 'rows');  % get unique rows
            d_star = u(:, 1)';
            for ff=1:nf  % deconstruct flows
                prop.(f{ff}) = u(:, 1 + ff)';
            end
            
            d = grid_i.edges{da_idx};  % points for integration
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_aac(d_star, d', prop, varargin{jj+1}{3:end});
            Lambda{ii} = permute(Lambda{ii}, [2,1]);

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, da_idx);
            [~,ll] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,ll,:);
            
            % [~,kk] = max(d_star == d_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);
            
            tools.textdone();
            
        
        %== BIN ==========================================================%
        %   When data input is binned (e.g., SP2 data).
        case {'bin', 'sp2'}
            disp(' Computing binned contribution ...');

            % Unpack inputs.
            s_idx = varargin{jj+1}{1};
            s_star2 = varargin{jj+1}{2};
            s_star = unique(s_star2)';
            s = grid_i.edges{s_idx};  % points for integration
            
            Lambda{ii} = full(tfer_bin(s_star', s'));
            
            % Duplicate over other grid dimensions.
            s2 = grid_i.elements(:, s_idx);
            [~,kk] = max(s == s2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            [~,kk] = max(s_star == s_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();

    end
end

% Loop over the various classifiers again to compile kernel.
disp(' Compiling kernel ...')
A = Lambda{1};  % initialize with first contribution
for ii=2:nc  % loop over other contributions
    A = A .* Lambda{ii};
end

% If second output selected, provide output prior to summing over charge.
if nargout > 1
    Ac = A;
end

A = sum(A, 3);  % sum over charge states

A = A .* grid_i.dr';  % multiply kernel by element area
A = sparse(A);  % exploit sparse structure in subsequent calculations

tools.textdone();
tools.textheader();

end


