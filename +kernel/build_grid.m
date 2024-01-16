
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

if mod(length(varargin), 2) ~= 0; error('Wrong number of inputs.'); end

tools.textheader('Computing kernel')

nc = (length(varargin)/2);  % number of classifiers
Lambda{nc} = [];  % initialize Lambda_ii

dm_idx = find(strcmp(grid_i.type, 'dm'));
mp_idx = find(strcmp(grid_i.type, 'mp'));
da_idx = find(strcmp(grid_i.type, 'da'));

if isempty(grid_i.type)  % set default indices, if type not specified (mass-mobility grid)
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
        prop_p = varargin{idx_p}{1};  % extract prop_pma from PMA input
        
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
        %   Computes charge fractions. The charger should be applied as the
        %   last contribution to ensure grid dimensions are respected. 
        case 'charger'
            disp(' Computing charger contribution...');
            
            % Assign inputs.
            if ~isempty(dm_idx)  % then PMA without DMA
                d = grid_i.edges{dm_idx};
                d2 = grid_i.elements(:, dm_idx);
                
            else  % otherwise inherit pre-computed value above
                d = dm';
                d2 = dm2';
            end
            
            f_z = charger(d', z_vec, varargin{jj+1}{1:end}); % get fraction charged for d vector
            Lambda{ii} = permute(f_z, [3, 2, 1]);
            
            % Duplicate over other grid dimensions.
            [~,kk] = max(d == d2, [], 2);
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
            if ~isempty(dm_idx)  % then PMA without DMA
                dm = grid_i.elements(:, dm_idx);
            end  % otherwise inherit pre-computed value above
            
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
            disp(' Computing AAC contribution...');
            
            addpath 'autils';
            
            % Assign inputs.
            d_star = grid_b.edges{da_idx};
            d = grid_i.edges{da_idx};  % points for integration
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_aac(d_star, d', varargin{jj+1}{1:end});
            Lambda{ii} = permute(Lambda{ii}, [2,1]);

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, da_idx);
            [~,kk] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            d_star2 = grid_b.elements(:, da_idx);
            [~,kk] = max(d_star == d_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);
            
            tools.textdone();
            
        
        %== BIN ==========================================================%
        %   When data input is binned (e.g., SP2 data).
        case {'bin', 'sp2'}
            disp(' Computing binned contribution...');

            % Unpack inputs.
            s_idx = varargin{jj+1}{1};
            s_star = grid_b.edges{s_idx};
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


