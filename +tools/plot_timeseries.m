
% PLOT_TIMESERIES  Plot timeseries data using a 2D histogram. 
%  Designed for CPMA-SP2 data. 
%  
%  IN is an input structure representing SP2 data, necessarily with fields
%  TimeDate and BHBL_BCmass. 
%  
%  AUTHOR: Timothy Sipkens, 2023-10-24

function [] = plot_timeseries(in, t_span, sp, left_cut)

if ~exist('sp', 'var'); sp = []; end
if ~exist('t_span', 'var'); t_span = []; end

if ~exist('left_cut', 'var'); left_cut = []; end
if isempty(left_cut); left_cut = 20; end

% Get unique minute-to-minute. 
ut = unique(in.TimeDate);

hx1 = logspace(log10(min(in.BHBL_BCmass)), ...
    log10(max(in.BHBL_BCmass)), 120);

% Enables cut of left edge of timeseries. 
idx_select = left_cut:length(ut);  % new indices with left cut

% Build histogram.
hy = zeros(length(idx_select),size(hx1,2)-1);
disp('Building histogram:')
tools.textbar([0,max(idx_select)]);
for ii=idx_select
    fl = in.TimeDate == ut(ii);
    hy(ii,:) = histcounts(in.BHBL_BCmass(fl), hx1);

    if mod(ii,10) == 0
        tools.textbar([ii,max(idx_select)]);
    end
end
tools.textbar([ii,max(idx_select)]);
hy(1:(idx_select(1)-1),:) = [];  % trim left cut points


% Plot results. 
batogram = hy';

figure(gcf);  clf;

% Plot histogram.
imagesc(idx_select, hx1(1:end-1), imgaussfilt(batogram, 1));

set(gca, 'YScale', 'log', 'YDir', 'normal');
ylim([2e-2, 3e1]);
ylabel('SP2 particle mass [fg]');
xlabel('Time index');

addpath('cmap'); colormap(tokyo);
colorbar;

% Plot overlay with selected data regions.
if ~isempty(sp)
    % Convert t_span from CPMA file to corresponding index. 
    idx_span1 = zeros(size(t_span,1),1);
    [tmp, ~] = find(ut == t_span(:,1)');  % convert t_span to equiv. idx
    idx_span1(any(ut == t_span(:,1)')) = tmp;
    
    idx_span2 = zeros(size(t_span,1),1);
    [tmp, ~] = find(ut == t_span(:,2)');  % convert t_span to equiv. idx
    idx_span2(any(ut == t_span(:,2)')) = tmp;

    if and(idx_span1(end) ~= 0, idx_span2(end) == 0)
        idx_span2(end) = length(ut);
    end
    
    hold on;
    xline(idx_span1, 'w:');
    xline(idx_span2, 'w:');
    hold off;
    
    hold on;
    plot(idx_span1, [sp.m_star] .* 1e18, 'ro');
    plot(idx_span2, [sp.m_star] .* 1e18, 'r>');
    plot([idx_span1, idx_span2]', [[sp.m_star]', [sp.m_star]']' .* 1e18, 'r-');
    hold off;
end
