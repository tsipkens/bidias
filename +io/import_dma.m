
% IMPORT_DMA  Import a raw CPC file containing SMPS data. 
% Author: Timothy Sipkens, 2020-04-12
%=========================================================================%

function [data,d_star,prop_dma,time_dma] = import_dma(fn)

prop_dma = kernel.prop_dma;

n = io.linecount(fn); % total line counts in file, sub-function below

opts = detectImportOptions(fn); % default options prior to modification below
opts.Whitespace = '\b ';
opts.Delimiter = {'\t'};



%-- Find header rows with specific data ----------------------------------%
optsc = opts;
optsc.DataLines = [1,(optsc.DataLines(1)+25)]; % adjust +25 in case detection failed 
optsc.VariableTypes(:) = {'string'};

ta = readtable(fn, optsc, 'ReadVariableNames', false);
tb = table2cell(ta);
n_date = find(cellfun(@(x) x=="Date",tb(:,1)));
n_r1 = find(cellfun(@(x) contains(x,"DMA Inner Radius"),tb(:,1)));
n_flow = find(cellfun(@(x) contains(x,"Sample Flow"),tb(:,1)));
n_T = find(cellfun(@(x) x=="Reference Gas Temperature (K)",tb(:,1)));
n_dia = find(cellfun(@(x) contains(x,"Diameter"),tb(:,1)));
%-------------------------------------------------------------------------%



tools.textheader('Reading SMPS file');  % output header indicating file processing

%-- Read header information ----------------------------------------------%
disp('Reading header information...');
disp('  Note: Info. stored in prop_dma should be checked.');
opts.VariableTypes(:) = {'double'};

opts.DataLines = [n_r1,n_r1+2];
ta = readtable(fn, opts, 'ReadVariableNames', false);
ta = table2array(ta);

% DMA inner radius (file is nominally normally in cm, covert to m).
prop_dma.R1 = ta(1,2);
if prop_dma.R1>0.1 % correct if illogical value for m, likely in cm then
    prop_dma.R1 = prop_dma.R1 / 100;
    disp('  Note: Converted R1 to m due to unexpected magnitude.');
end
disp(['    R1 = ', num2str(prop_dma.R1 * 100),' cm']);


% DMA outer radius.
prop_dma.R2 = ta(2,2);
if prop_dma.R2>0.1 % correct if illogical value for m, likely in cm then
    prop_dma.R2 = prop_dma.R2 / 100;
    disp('   Note: Converted R2 to m due to unexpected magnitude.');
end
disp(['    R2 = ', num2str(prop_dma.R2 * 100),' cm']);


% DMA column length. 
prop_dma.L = ta(3,2);
if prop_dma.L>1 % correct if illogical value for m, likely in cm then
    prop_dma.L = prop_dma.L / 100;
    disp('  Note: Converted L to m due to unexpected magnitude.');
end
disp(['    L = ', num2str(prop_dma.L * 100),' cm']);


% Read flow rates.
if ~isempty(n_flow)
    opts.DataLines = [n_flow,n_flow]; % this is unreliable (i.e. row may vary)
    ta = readtable(fn, opts, 'ReadVariableNames', false);
    ta = table2array(ta);
    prop_dma.Q_a = ta(1,2)/60/1000;
    prop_dma.Q_s = prop_dma.Q_a; % equal flow assumption
    prop_dma.Q_c = ta(1,4)/60/1000;
    prop_dma.Q_m = prop_dma.Q_c; % equal flow assumption
    disp('  Note: Applied equal flow assumption for Q_s and Q_m.');
    disp('  Note: Reading flow information from the header is unreliable.');
end
disp(['    Q_a = ', num2str(prop_dma.Q_a),' m3/s = ', ...
    num2str(prop_dma.Q_a*60*1000), ' LPM']);
disp(['    Q_c = ', num2str(prop_dma.Q_c),' m3/s = ', ...
    num2str(prop_dma.Q_c*60*1000), ' LPM']);


% Read temperature and pressure.
if ~isempty(n_T)
    opts.DataLines = [n_T,n_T+1];
    ta = readtable(fn, opts, 'ReadVariableNames', false);
    ta = table2array(ta);
    prop_dma.T = ta(1,2); % more DMA properties
    prop_dma.p = ta(2,2)/101.325;
    
    disp(['    T = ', num2str(prop_dma.T),' K']);
    disp(['    p = ', num2str(prop_dma.p),' atm']);
end
disp('Complete.');
disp(' ');
%-------------------------------------------------------------------------%



%-- Get number of samples/data columns -----------------------------------%
opts.DataLines = [40, 40]; % row 40 chosen arbitrarily to avoid header
ta = readtable(fn, opts, 'ReadVariableNames', false);

t0 = cellfun(@(x) isnan(x),table2cell(ta(1,:)),'UniformOutput',false);
idx_endc = find([t0{:}],1); % remove trailing columns
if isempty(idx_endc); idx_endc = width(ta); end

ta = table2array(ta);
ncol = idx_endc-1; % number of data samples/scans
%-------------------------------------------------------------------------%



%-- Read date/time rows --------------------------------------------------%
disp('Reading date...');
optsb = opts;
optsb.VariableTypes(:) = {'string'};
optsb.DataLines = [n_date,n_date+1];

ta = readtable(fn, optsb, 'ReadVariableNames', false);
ta = table2array(ta(:,1:idx_endc));
for ii=1:ncol
    time_dma(ii,1) = datetime(ta(1,ii+1),'InputFormat','MM/dd/yy')+...
        duration(ta(2,ii+1),'InputFormat','hh:mm:ss');
end
disp('Complete.');
disp(' ');
%-------------------------------------------------------------------------%



%-- Read in actual data --------------------------%
disp('Reading data...');
idx_0 = n_dia; % in case automatic detection started at date line

for ii=0:10 % loop through several rows to find a number
    idx_0 = idx_0 + 1;
    opts.DataLines = [idx_0,idx_0];
    ta = readtable(fn, opts, 'ReadVariableNames', false);
    if ~isnan(table2array(ta(1,1))); break; end
end

opts.DataLines = [idx_0,n-26]; % read in the data
ta = readtable(fn, opts, 'ReadVariableNames', false);
idx_0 = table2array(ta(:,1:idx_endc));
data = idx_0(:,2:(ncol+1));
d_star = idx_0(:,1);
%-------------------------------------------------%

tools.textheader();


%-- Set remaining prop_dma parameters --------%
%   Note: update these after this function call, if necessary.
prop_dma.G_DMA = 3.0256;
prop_dma.bet = 0.1000;
prop_dma.del = 0;
%---------------------------------------------%


end

