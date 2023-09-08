baseline_tag = {'noodor', 'baseline', 'spont.'};
blank_tag = {'ACSF'};
this_group = {'trained1','trained2'};
cat.traces = [];
todo_fish = [];
for i_fish = 1:numel(experiment.series)
    if ismember(experiment.series{i_fish}.group,this_group)
        todo_fish = [todo_fish;i_fish];
    end
end

% Initialize odor-trials indices
reference_data = experiment.series{todo_fish(1)}.data;
is_odor_trial = ~ ( ismember(reference_data.stim_type, baseline_tag) | ...
                    ismember(reference_data.stim_type, blank_tag));

% ... and sparseness matrix [trials x fish]
population_sparsenessMAT = nan(sum(is_odor_trial), ...
                               numel(todo_fish));

%% Calculate and store sparseness for each trial and fish

for i_fish = 1:numel(todo_fish)
    data = experiment.series{todo_fish(i_fish)}.data;
    fs = data.meta.framerate;
    interval = [floor(data.stim_on_sec*fs) : 1+floor(data.stim_off_sec*fs)];

    % extract activity traces (only for odor periods)
    activityTraces = traceFormat(selectTimeFrame(data.tracesdn,interval,data.L),length(interval));
    
    % eliminate any baseline or blank trials
    activityTraces(:,:,~is_odor_trial) = [];
    
    % calculate sparseness
    population_sparseness = calculateSparseness(activityTraces,'population');

    % store sparseness in output matrix
    population_sparsenessMAT(:,i_fish) = population_sparseness;
end

% data for plotting
dataMAT = population_sparsenessMAT;
labs = reference_data.stim_type(is_odor_trial);

%% plot over all odor trials (in chronological order)

data = dataMAT;

t = 1:size(data,1);
y = nanmean(data,2);
figure; hold on
curve1 = y + std(data,[],2,'omitnan'); curve2 = y - std(data,[],2,'omitnan');
h = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'r','FaceAlpha',.3,'EdgeColor','none')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y, 'LineWidth', 5, 'Color', 'r')

xticks(t); xticklabels(labs)
xlabel('')
ylabel('pop. sparseness')
axis tight
ylim([0, 1])

set(gcf, 'Position', [50 50 400 280]);
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 


%% plot over all odor trials (in stimulus order)

[~,idx] = sort(labs);
data = dataMAT(idx,:);

t = 1:size(data,1);
y = nanmean(data,2);
figure; hold on
curve1 = y + std(data,[],2,'omitnan'); curve2 = y - std(data,[],2,'omitnan');
h = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'r','FaceAlpha',.3,'EdgeColor','none')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y, 'LineWidth', 5, 'Color', 'r')

xticks(t); xticklabels(labs(idx))
xlabel('')
ylabel('pop. sparseness')
axis tight
ylim([0, 1])

set(gcf, 'Position', [50 50 400 280]);
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 

%% plot over relative trials of each odor

mode = 'all';

switch mode
    case 'all'
        stims = unique(labs);
    case 'trained only'
        stims = {'Ala','Arg','His'};
    case 'naive only'
        stims = {'Leu','Ser','Trp'};
    otherwise
        stims = {};
end



% find the maximum number of repetions for any stimulus
maxreps = 0;
for i_stim = 1:numel(stims)
    thisstim = stims(i_stim);
    idx = ismember(labs,thisstim);
    maxreps = max([maxreps, sum(idx)]);
end

% define data: trials x (stims*fish)
data = nan(maxreps,numel(stims)*numel(todo_fish));
for i_stim = 1:numel(stims)
    thisstim = stims(i_stim);
    idx = ismember(labs,thisstim);
    bias = (i_stim-1)*numel(todo_fish);
    data(1:sum(idx),bias+[1:numel(todo_fish)]) = dataMAT(idx,:);
end

% plot
t = 1:size(data,1);
y = nanmean(data,2);
figure; hold on
curve1 = y + std(data,[],2,'omitnan'); curve2 = y - std(data,[],2,'omitnan');
h = patch([t,fliplr(t)],[curve1; fliplr(curve2')'],'r','FaceAlpha',.3,'EdgeColor','none')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t,y, 'LineWidth', 5, 'Color', 'r')

xticks(t);
xlabel('')
ylabel('pop. sparseness')
axis tight
ylim([0, 1])

set(gcf, 'Position', [50 50 400 280]);
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 

% stats
ranksum(data(1,:),data(end,:))


