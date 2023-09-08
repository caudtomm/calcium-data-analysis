function [avg_d, sem_d, std_d, C] = plotDistMatInTime(data,interval_len,interval_step, cells, trialn2use, pl, s, cax, tag, method)
% Distance between trials on a sliding time window basis (using avg 
% activity vector)
% 
% [avgd, semd, stdd] = plotCosDMatInTime(data,interval_len,interval_step,
%   cells,trialn2use,pl,s,cax,tag,method)
%
% ex.
% data = data; (struct from _DATA file (see ProcessingPipeline))
% interval_len = 5; % in sec
% interval_step = 3; % in sec
% cells = [1;2;3;4]; % [N:1] cell IDs
% trialn2use = [1;4;5]; % [n:1] trials to do
% pl = true; : plot figures, boolean
% s = true; : save figures, boolean
% cax = [.7 1]; % caxis of plots
% tag = 'trialnum'; % attached to autput file name
% method = 'cosine'; % which distance metric to use
%
%


% initialize vars
grouplab = data.meta.seriesid;
traces = data.tracesdn(:,trialn2use);
L = data.L;
traces = selectCells(traces, L, cells');
fs = data.meta.framerate;
interval = 1:(interval_len*fs+1);

avg_d = [];
sem_d = [];
std_d = [];
C = [];

% create output directory
PathOut = fullfile(grouplab, 'figures', mfilename,'tmp');
if s && ~exist(PathOut, 'dir'); mkdir(PathOut); end

% loop through time bins
for i_t = 1:floor((L-length(interval))/(interval_step*fs))
    
    % specify frames belonging to the current time bin
    thisinterval = interval+floor((i_t-1)*interval_step*fs);
    
    % compute distance matrix
    c = getDistinRegn(traces,L,thisinterval,method);
    
    if pl
        % plotting
        figure; imagesc(1-c); axis square; hold on
        xticks(1:size(c,1)); xticklabels(data.stim_type(idx)); xtickangle(90)
        yticks(1:size(c,1)); yticklabels(data.stim_type(idx))
        ylabel('Stimulus type')
        title(['1 - ',method,': sec ', ...
            num2str(thisinterval(1)/fs-data.stim_on_sec), '-', ...
            num2str(thisinterval(end)/fs-data.stim_on_sec)])
        clim(cax); colorbar
        hold off

        % saving figure
        if s; saveas(gcf, [PathOut, '\', method, '_', tag, '_', ...
                num2str(i_t,'%04d'), '.tif'], 'tif'); end
    end
    
    % store and return vars
    c = c(triu(c)~=0);
    avg_d = [avg_d; nanmean(c)];
    sem_d = [sem_d; std(c)/sqrt(numel(c))];
    std_d = [std_d; std(c)];
    % C : [PxB], where:
    % P is the num of elements in the upper triangle of the 
    % single-time-bin intertrial distance matrix (without diagonal); 
    % B is the number of time bins
    C(:,i_t) = c;
    
    
end
    
    
    
    
    