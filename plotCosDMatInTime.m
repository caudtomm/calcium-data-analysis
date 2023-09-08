function [avg_cosd, sem_cosd, std_cosd, C] = plotDMatInTime(data,interval_len,interval_step, cells, idx, s, cax, tag, method)
% Cosine distance between trials on a sliding time window basis (using avg 
% activity vector)
% 
% plotCosDMatInTime(data,interval_len,interval_step, s)
% 
% interval_len and interval_step in seconds
%
% ex.
% data = data; (struct from _DATA file (see ProcessingPipeline))
% interval_len = 5; % in sec
% interval_step = 3; % in sec
% cells = [1:10]'; % can be extracted from data.anat_regions.cells
% s = true; : save figures, boolean
% cax = [.7 1]; % caxis of plots
%


% grouplab = 'TC_200330_TC0004_01_sxpDp';
grouplab = data.meta.seriesid;

% traces = data.tracesdns(:,idx);
% L = data.Ldns;

traces = data.tracesdn(:,idx);
L = data.L;

traces = selectCells(traces, L, cells');
n = numel(cells)

% fs = data.meta.framerate / data.meta.downsample;
fs = data.meta.framerate;

interval = 1:(interval_len*fs+1);

avg_cosd = [];
sem_cosd = [];
std_cosd = [];
C = [];

PathOut = [grouplab, '\figures\', mfilename,'\tmp'];
if s && ~exist(PathOut, 'dir'); mkdir(PathOut); end

for i_t = 1:floor((L-length(interval))/(interval_step*fs))
    
    thisinterval = interval+floor((i_t-1)*interval_step*fs);
    
    c = getDistinRegn(traces,L,thisinterval,method);
    
    
    figure; imagesc(1-c); axis square; hold on
    xticks(1:size(c,1)); xticklabels(data.stim_type(idx)); xtickangle(90)
    yticks(1:size(c,1)); yticklabels(data.stim_type(idx))
    xlabel('Trial number (ITI: 180 s)')
    ylabel('Stimulus type')
    title(['1 - ',method,': sec ', ...
        num2str(thisinterval(1)/fs-data.stim_on_sec), '-', ...
        num2str(thisinterval(end)/fs-data.stim_on_sec)])
    caxis(cax); colorbar
    hold off
    
    if s; saveas(gcf, [PathOut, '\', method, '_', tag, '_', ...
            num2str(i_t,'%04d'), '.tif'], 'tif'); end
    
    c = c(triu(c)~=0);
    avg_cosd = [avg_cosd; nanmean(c)];
    sem_cosd = [sem_cosd; std(c)/sqrt(numel(c))];
    std_cosd = [std_cosd; std(c)];
    C(:,i_t) = c;
    
    
end
    
    
    
    
    