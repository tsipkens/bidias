
% GEN_GENERAL  A general function to compute a 2D transfer function. 
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

function A = gen_general(grid_i, z_vec, varargin)

tools.textheader('Computing kernel')

nc = (length(varargin)/2);  % number of classifiers
Lambda{nc} = [];  % initialize Lambda_ii

dm_idx = find(strcmp(grid_i.type, 'dm'));
mp_idx = find(strcmp(grid_i.type, 'mp'));

if isempty(dm_idx)  % set default indices, if type not specified
    dm_idx = 2;
    mp_idx = 1;
end

% Loop over the various classifiers. 
for ii=1:nc
    jj = 2*ii - 1;  % index of classifier

    switch varargin{jj}
        
        %== CHARGER ======================================================%
        case 'charger'
            disp(' Computing charger contribution...');
            
            % Assign inputs.
            d = grid_i.edges{dm_idx};
            f_z = charger(d', z_vec, varargin{jj+1}{1:end}); % get fraction charged for d vector
            
            Lambda{ii} = permute(f_z, [3, 2, 1]);

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, dm_idx);
            [~,kk] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            tools.textdone();


        %== SMPS =========================================================%
        case {'dma', 'smps'}
            disp(' Computing DMA contribution...');
            
            % Unpack inputs.
            d_star2 = varargin{jj+1}{1};
            d_star = unique(d_star2)';
            d = grid_i.edges{dm_idx};  % points for integration
            prop_dma = varargin{jj+1}{2};  % DMA properties
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_dma(d_star, d', z_vec, prop_dma);

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
            disp(' Computing PMA contribution...');

            % Unpack inputs.
            % Points for integration.
            r = grid_i.elements;  % get elements from the grid
            m = r(:, mp_idx);  % masses at which to compute the transfer function (not setpoints)
            d = r(:, dm_idx);  % mobilities at which to compute the transfer function (not setpoints)
            
            % Unpack inputs.
            m_star = varargin{jj+1}{1};  % don't use unique(), as resolution may change
            prop_p = varargin{jj+1}{2};  % DMA properties

            addpath 'tfer\tfer-pma';  % added to calculate sp
            sp = get_setpoint(prop_p,...  % get PMA setpoints
                'm_star', m_star .* 1e-18, ...  % mass from the grid
                varargin{jj+1}{3:end});  % extra name-value pair to specify setpoint
            
            Lambda{ii} = tfer_pma(...
                sp, m, d, z_vec, prop_p);

            tools.textdone();
            

        %== AAC ==========================================================%
        case 'aac'
            
        
        %== BIN ==========================================================%
        %   When data input is binned (e.g., SP2 data).
        case 'bin'
            disp(' Computing binned contribution...');

            % Unpack inputs.
            Lambda{ii} = 1;

            % Duplicate over other grid dimensions.
            m_star2 = grid_b.elements(:, ii);
            [~,kk] = max(m_star == m_star2, [], 2);
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


