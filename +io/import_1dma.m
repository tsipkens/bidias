
% IMPORT_1DMA  Import a raw CPC file contianing SMPS data. 
% Author: Timothy Sipkens, 2020-04-12
%=========================================================================%

function [data,d_star,prop_dma] = import_1dma(fn)

n = io.linecount(fn);
opts = detectImportOptions(fn);
opts.Whitespace = '\b ';
opts.Delimiter = {'\t'};

optsc = opts;
idx_0 = optsc.DataLines(1);
for ii=1:idx_0
    optsc.DataLines = [ii,ii];
    for jj=1:length(optsc.VariableTypes)
        optsc.VariableTypes{jj} = 'string';
    end
    
    ta = readtable(fn,optsc);
    if or(table2array(ta(1,1))=="DMA Inner Radius(cm)",table2array(ta(1,1))=="DMA Inner Radius (cm)")
        nR1 = ii;
    end
    if table2array(ta(1,1))=="Reference Gas Temperature (K)"; nT = ii; end
    if table2array(ta(1,1))=="Date"; nDate = ii; end
end

opts.DataLines = [40,40]; % use to get number of samples
opts.VariableTypes(:) = {'double'};
ta = readtable(fn,opts);

% idx_exclude = isnan(table2array(ta(:,1)));
% idx_end = find(idx_exclude,1); % remove trailing rows

t0 = cellfun(@(x) isnan(x),table2cell(ta(1,:)),'UniformOutput',false);
idx_endc = find([t0{:}],1); % remove trailing columns

% ta = ta(1:(idx_end-1),1:(idx_endc-1));

ta = table2array(ta);
ncol = length(ta(:,1:(idx_endc-1)))-1; % number of data samples/scans

opts.DataLines = [nR1,nR1+2];
ta = readtable(fn,opts);
ta = table2array(ta(:,1:(idx_endc-1)));
prop_dma.R1 = ta(1,2);
prop_dma.R2 = ta(2,2);
prop_dma.L = ta(3,2); % some DMA properties

opts.DataLines = [12,12]; % this is unreliable
ta = readtable(fn,opts);
ta = table2array(ta(:,1:(idx_endc-1)));
prop_dma.Q_a = ta(1,2)/60/1000;
prop_dma.Q_s = prop_dma.Q_a; % equal flow assumption
prop_dma.Q_c = ta(1,4)/60/1000;
prop_dma.Q_m = prop_dma.Q_c; % equal flow assumption

opts.DataLines = [nT,nT+1];
ta = readtable(fn,opts);
ta = table2array(ta(:,1:(idx_endc-1)));
prop_dma.T = ta(1,2); % more DMA properties
prop_dma.p = ta(2,2);


optsb = opts;
for ii=1:length(optsb.VariableTypes); optsb.VariableTypes{ii} = 'string'; end
optsb.DataLines = [nDate,nDate+1];
ta = readtable(fn,optsb);
ta = table2array(ta(:,1:(idx_endc-1)));
for ii=1:ncol
    date_smps(ii,1) = datetime(ta(1,ii+1),'InputFormat','MM/dd/yy')+...
        duration(ta(2,ii+1),'InputFormat','hh:mm:ss');
end


opts.DataLines = [idx_0,n-26]; % read in the data
ta = readtable(fn,opts);
idx_0 = table2array(ta(:,1:(idx_endc-1)));
data = idx_0(:,2:(ncol+1));
d_star = idx_0(:,1);

end

