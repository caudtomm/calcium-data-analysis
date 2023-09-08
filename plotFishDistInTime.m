function [c,interval] = plotFishDistInTime(data,secrange,method,s)

%data = experiment.series{5}.data
traces = data.tracesdn;
traces = traceFormat(traces,data.L);
ntrials = numel(data.trial_num);

if ~exist('secrange','var')
    secrange = [-5, 15]; % from [stim_on, stim_off]
end
if ~exist('method','var')
    method = 'cosine'
end


fs = data.meta.framerate;
interval =  floor((data.stim_on_sec+secrange(1))*fs) : ...
            floor((data.stim_off_sec+secrange(2))*fs);

for i_plot = 1:2
    h = figure;
    
    switch i_plot
        case 1
            ttl = 'by trialnum';
            [~,idx] = sort(data.trial_num);
        case 2
            ttl = 'by stimtype';
            idx = data.idx_by_stim_type;
        otherwise
            error('plotting function out of range')
    end

    % computing
    traces2 = nan(data.L*ntrials,data.N);
    for i=1:ntrials
    traces2(1+data.L*(i-1):data.L*i,:) = traces(:,:,idx(i));
    end
    traces2 = selectTimeFrame(traces2,interval,data.L);
    c= squareform(pdist(traces2,method));

    % plotting
    t = [0:1/fs:size(c,1)/fs-1/fs]';
    T = data.stim_off_sec - data.stim_on_sec + diff(secrange);
    imagesc(t/T,t/T,1-c)

    % cosmetics and labels
    axis square
    col = colorbar('Color','k');
    ylabel(col,'cosine similarity','Color','k')
    xticks(0:ntrials-1); xticklabels(data.stim_type(idx)); xtickangle(90)
    yticks(0:ntrials-1); yticklabels(data.stim_type(idx));
    ylabel('trials','Color','k')
    xlabel('trials','Color','k')
    caxis([.2 1])
    colormap(['parula'])
    set(gca, 'color', 'none', 'XColor','k', 'YColor','k', 'ZColor','k');
    set(gcf, 'color', 'w'); 
    title(ttl,'Color','k')

    % saving
    pathout = fullfile('figures');
    FileOut = fullfile(pathout,strcat(mfilename,ttl,'_',method,'.fig'));
    if s; savefig(h,FileOut); end
end