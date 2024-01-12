
% BUILD_GRID  A general function to compute a 2D transfer function. 
%  
%  NOTE: This does not apply to PMA-DMA systems, where mobility diameter 
%  has an impact on the PMA transfer function. 
%  
%  Input is name value pairs for classifier/transfer function type and
%  a cell of inputs corresponding to the transfer function evaluation. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2024-01-11

function A = build_grid(grid_b, grid_i, z_vec, varargin)

tools.textheader('Computing kernel')

nc = (length(varargin)/2);  % number of classifiers
Lambda{nc} = [];  % initialize Lambda_ii

dm_idx = find(strcmp(grid_i.type, 'dm'));
mp_idx = find(strcmp(grid_i.type, 'mp'));

if isempty(grid_i.type)  % set default indices, if type not specified (mass-mobility grid)
    dm_idx = 2;
    mp_idx = 1;
end

% Loop over the various classifiers. 
for ii=1:nc
    jj = 2*ii - 1;  % index of classifier

    switch varargin{jj}
        
        %== CHARGER ======================================================%
        %   Computes charge fractions. The charger should be applied as the
        %   last contribution to ensure grid dimensions are respected. 
        case 'charger'
            disp(' Computing charger contribution...');
            
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
        %   NOTE: The DMA transfer function is 1D (only a function of 
        %   mobility), which is exploited to speed evaluation. 
        case {'dma', 'smps'}
            disp(' Computing DMA contribution...');
            
            % Assign inputs.
            d_star = grid_b.edges{dm_idx};  % DMA setpoints
            d = grid_i.edges{dm_idx};  % points for integration
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_dma(d_star, d', z_vec, varargin{jj+1}{1:end});

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, dm_idx);
            [~,kk] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            d_star2 = grid_b.elements(:, dm_idx);
            [~,kk] = max(d_star == d_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();

        
        %== PMA ==========================================================%
        %   Currently assumes other dimension is mobility diameter. 
        case 'pma'
            disp(' Computing PMA contribution...');

            % Unpack inputs.
            m_star = grid_b.edges{mp_idx};  % DMA setpoints
            prop_p = varargin{jj+1}{1};  % DMA properties

            % Points for integration.
            m = grid_i.elements(:, mp_idx);  % masses at which to compute the transfer function (not setpoints)
            if isempty(dm_idx)  % then PMA without DMA
                dm = (m .* 1e-18 ./ prop_p.rho0) .^ ...
                    (1./prop_p.Dm) .* 1e9;  % use mass-mobility
                disp('  Invoking mass-mobility relation for PMA.')
            else
                dm = grid_i.elements(:, dm_idx);
            end
            
            addpath 'tfer\tfer-pma';  % added to calculate sp
            sp = get_setpoint(prop_p,...  % get PMA setpoints
                'm_star', m_star .* 1e-18, ...  % mass from the grid
                varargin{jj+1}{2:end});  % extra name-value pair to specify setpoint
            
            Lambda{ii} = tfer_pma(...
                sp, m, dm, z_vec, prop_p);
            
            % Duplicate over other grid dimensions.
            m_star2 = grid_b.elements(:, mp_idx);
            [~,kk] = max(m_star == m_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);
    
            tools.textdone();
            

        %== AAC ==========================================================%
        case 'aac'
            
        
        %== BIN ==========================================================%
        %   When data input is binned (e.g., SP2 data).
        case {'bin', 'sp2'}
            disp(' Computing binned contribution...');

            % Unpack inputs.
            s_idx = varargin{jj+1}{1};  % DMA properties
            s_star = grid_b.edges{s_idx};  % DMA setpoints
            s = grid_i.edges{s_idx};  % points for integration
            
            Lambda{ii} = full(tfer_bin(s_star', s'));
            
            % Duplicate over other grid dimensions.
            s2 = grid_i.elements(:, s_idx);
            [~,kk] = max(s == s2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            s_star2 = grid_b.elements(:, s_idx);
            [~,kk] = max(s_star == s_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();

    end
end

% Loop over the various classifiers again to compile kernel.
A = Lambda{1};  % initialize with first contribution
for ii=2:nc  % loop over other contributions
    A = A .* Lambda{ii};
end
A = sum(A, 3);  % sum over charge states

A = A .* grid_i.dr';  % multiply kernel by element area
A = sparse(A);  % exploit sparse structure in subsequent calculations

tools.textheader();

end


