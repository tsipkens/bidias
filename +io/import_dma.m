
% IMPORT_DMA  Import a raw CPC file contianing SMPS data. 
% Author: Timothy Sipkens, 2020-04-12
%=========================================================================%

function [data,d_star,prop_dma,time_dma] = import_dma(fn)

n = io.linecount(fn); % total line counts in file, sub-function below

opts = detectImportOptions(fn); % default options prior to modification below
opts.Whitespace = '\b ';
opts.Delimiter = {'\t'};


%-- Find header rows with specific data ----------------------------------%
optsc = opts;
optsc.DataLines = [1,(optsc.DataLines(1)+9)];
optsc.VariableTypes(:) = {'string'};

ta = readtable(fn,optsc);
tb = table2cell(ta);
n_date = find(cellfun(@(x) x=="Date",tb(:,1)));
n_r1 = find(or(...
    cellfun(@(x) x=="DMA Inner Radius(cm)",tb(:,1)),...
    cellfun(@(x) x=="DMA Inner Radius (cm)",tb(:,1))));
n_T = find(cellfun(@(x) x=="Reference Gas Temperature (K)",tb(:,1)));

idx_0 = n_date + 8; % in case automatic detection started at date line
%-------------------------------------------------------------------------%


%-- Read header information ----------------------------------------------%
opts.VariableTypes(:) = {'double'};

opts.DataLines = [n_r1,n_r1+2];
ta = readtable(fn,opts);
ta = table2array(ta);
prop_dma.R1 = ta(1,2)/100; % some DMA properties (file is in cm, covert to m)
prop_dma.R2 = ta(2,2)/100;
prop_dma.L = ta(3,2)/100;

opts.DataLines = [12,12]; % this is unreliable (i.e. row may vary)
ta = readtable(fn,opts);
ta = table2array(ta);
prop_dma.Q_a = ta(1,2)/60/1000;
prop_dma.Q_s = prop_dma.Q_a; % equal flow assumption
prop_dma.Q_c = ta(1,4)/60/1000;
prop_dma.Q_m = prop_dma.Q_c; % equal flow assumption

opts.DataLines = [n_T,n_T+1];
ta = readtable(fn,opts);
ta = table2array(ta);
prop_dma.T = ta(1,2); % more DMA properties
prop_dma.p = ta(2,2)/101.325;
%-------------------------------------------------------------------------%


%-- Get number of samples/data columns -----------------------------------%
opts.DataLines = [40,40]; % row 40 chosen arbitrarily to avoid header
ta = readtable(fn,opts);

t0 = cellfun(@(x) isnan(x),table2cell(ta(1,:)),'UniformOutput',false);
idx_endc = find([t0{:}],1); % remove trailing columns
if isempty(idx_endc); idx_endc = width(ta); end

ta = table2array(ta);
ncol = idx_endc-1; % number of data samples/scans
%-------------------------------------------------------------------------%


%-- Read date/time rows --------------------------------------------------%
optsb = opts;
optsb.VariableTypes(:) = {'string'};
optsb.DataLines = [n_date,n_date+1];

ta = readtable(fn,optsb);
ta = table2array(ta(:,1:idx_endc));
for ii=1:ncol
    time_smps(ii,1) = datetime(ta(1,ii+1),'InputFormat','MM/dd/yy')+...
        duration(ta(2,ii+1),'InputFormat','hh:mm:ss');
end
%-------------------------------------------------------------------------%


%-- Read in actual data --------------------------%
opts.DataLines = [idx_0,n-26]; % read in the data
ta = readtable(fn,opts);
idx_0 = table2array(ta(:,1:idx_endc));
data = idx_0(:,2:(ncol+1));
d_star = idx_0(:,1);
%-------------------------------------------------%


%-- Set remaining prop_dma parameters --------%
%   Note: update these after this function call, if necessary.
prop_dma.G_DMA = 3.0256;
prop_dma.bet = 0.1000;
prop_dma.del = 0;
%---------------------------------------------%


end

