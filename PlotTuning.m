% extract tuning curves
% for all cells in 'cells'
% activity and avg across stimulus types in data.stim_type
% avg during stim window


L = data.Ldns;
traces = selectCells(data.tracesdns,L,cells');
fs = data.meta.framerate/data.meta.downsample;
x = 0:1/fs:L/fs-1/fs;
% interval = floor(data.stim_on_sec*fs):floor(data.stim_off_sec*fs);
interval = 1:10;


% get only relevant time window
traces = selectTimeFrame(traceFormat(traces,L),interval,L);
% now traces : [time x cells x trials]

% get average activity per cell, per trial
tuning = squeeze(nanmean(traces,1));
% tuning : [cells x trials]

% get avg tuning per stim type
allstimtypes = unique(data.stim_type);
tuning_avg = zeros(size(tuning,1),numel(allstimtypes));
for i_stim = 1:numel(allstimtypes)
    
    thisstim = allstimtypes{i_stim};
    thisstimids = ismember(data.stim_type,thisstim);
    
    %get only trials of this type
    tmpmat = tuning(:,thisstimids);
    
    % average and store
    tuning_avg(:,i_stim) = nanmean(tmpmat,2);
    
end


%% plotting 

% plot average tuning
figure;
imagesc(tuning_avg);
title('average responses')
ylabel('cell #')
xticklabels(allstimtypes)


% plot tuning per cell
for i_cell = 1:10%numel(cells)

    i_cell = randi(numel(cells),1);
    
    figure; hold on

    % plot single-trial activity
    for i_stim = 1:numel(allstimtypes)

        thisstim = allstimtypes{i_stim};
        thisstimids = ismember(data.stim_type,thisstim);

        %plots
%             boxplot(tuning(i_cell,thisstimids))
        scatter(i_stim*ones(1,sum(thisstimids)) ...
                + .1*randn(1,sum(thisstimids)), ...
                tuning(i_cell,thisstimids))

    end

    % plot average activity per stimulus
    scatter(1:numel(allstimtypes), ...
            tuning_avg(i_cell,:),50,'filled')


   % cosmetics
   title(['cell #', num2str(cells(i_cell))])
   xticks(1:numel(allstimtypes))
   xticklabels(allstimtypes)
   xlim([.3, numel(allstimtypes)+.7])
   ylabel('average dF/F')

end


% manhattan plot
figure; hold on
for i_stim = 1:numel(allstimtypes)

    scatter(i_stim*ones(1,numel(cells)), ...
        tuning_avg(:,i_stim),10,'k')

end
% superimposed violin plot
violin(tuning_avg)
% cosmetics
title('average responses across trials')
xticks(1:numel(allstimtypes))
xticklabels(allstimtypes)
xlim([.3, numel(allstimtypes)+.7])
ylabel('average dF/F')
    
    



