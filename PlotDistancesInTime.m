
%% DISTANCE IN TIME BETWEEN TRIALS OF THE THE SAME VS DIFFERENT ODORS

todo_fish = 1:numel(experiment.series);

stims{1} = {'spont.', 'baseline', 'noodor'};
stims{2} = {'ACSFodor','ACSF'};
stims{3} = {'Hodor','Trp'};
stims{4} = {'Sodor','Ser'};
stims{5} = {'Aodor','Ala'};
stims{6} = {'FOodor','Food'};
%stims{7} = {'SXodor'};

stims2use = [3:6];
cax = [.60 1];
method = 'cosine';

% time interval of interest, in seconds:
% [from stim_on, from stim_off]
rel_time = [-5, 19];

labs = {'avg diff. odor', 'avg same odor'};

allDsamestim = [];
allDdiffstim = [];
for i_fish = 1:numel(todo_fish)
    data = experiment.series{todo_fish(i_fish)}.data;
    
    % select cells
    cells = [1:data.N]';
    
    % select trials
    idx = [];
    for i_stim = 1:numel(stims2use)
        thisstim = stims2use(i_stim);
        idx2 = find(ismember(data.stim_type,stims{thisstim}));
        idx = [idx; idx2'];
    end
    trialn2use = sort(idx);

    %% get all intertrial distances
    
    % get intertrial distance matrix FOR THE WHOLE TRIAL
    % C : [PxB], where:
    % P is the num of elements in the upper triangle of the 
    % single-time-bin intertrial distance matrix (without diagonal); 
    % B is the number of time bins
    [~,~,~,C] = plotDistMatInTime(data,2,1,cells,trialn2use,0,0,cax,'trialnum',method);
    
    % define the avg distance in the first 4 bins to be equal to zero (we
    % are measuring DELTA distance)
    fs = data.meta.framerate;
    interval = [round(data.stim_on_sec)+rel_time(1) : ...
        min([round(data.stim_off_sec)+rel_time(2), size(C,2)])]';
    %C = C - nanmean(C(:,interval(1)+[0:3]),2);

    % identify trials of the same stimulus
    idx = [];
    for i_stim = 1:numel(stims2use)
        thisstim = stims2use(i_stim);
        idx2 = find(ismember(data.stim_type(trialn2use),stims{i_stim}));
        idx = [idx, getLinearCorrespondanceMap(idx2)];
    end
    idx = sort(idx);
    issamestim = zeros(size(C,1),1);
    issamestim(idx) = 1; issamestim = logical(issamestim);

    % store the fish's average in the time interval of interest
    allDsamestim = [allDsamestim;nanmean(C(issamestim,interval),1)]; % distances btw trials of the same stimulus
    allDdiffstim = [allDdiffstim;nanmean(C(~issamestim,interval),1)]; % distances btw trials of different stimuli
end

%% FIGURE: distances across trials in individual fish
cscale = .3;
figure
subplot(131); imagesc(allDsamestim); title('dist same stim'); ylabel(['# fish']);% caxis(cscale*[-1,1]);
subplot(132); imagesc(allDdiffstim); title('dist diff stim');% caxis(cscale*[-1,1])
subplot(133); plot(nanmean(allDdiffstim)); hold on; plot(nanmean(allDsamestim));axis tight;title('avg'), legend(labs);%ylim(cscale*[-1,1])
    
%% FIGURE: distances across trials for all fish
t = linspace(-5,39,length(interval));
figure;

y0 = allDdiffstim;
y = nanmean(y0,1);
curve1 = y + std(y0,[],1,'omitnan');
curve2 = y - std(y0,[],1,'omitnan');
patch([t,fliplr(t)],[curve1, fliplr(curve2)],'k','EdgeColor','none','FaceAlpha',.3)
hold on
b(1) = plot(t,y,'Color','k','LineWidth',3);

%plot(t,y0,'k')

y0 = allDsamestim;
y = nanmean(y0,1);
curve1 = y + std(y0,[],1,'omitnan');
curve2 = y - std(y0,[],1,'omitnan');
patch([t,fliplr(t)],[curve1, fliplr(curve2)],'r','EdgeColor','none','FaceAlpha',.3)
b(2) = plot(t,y,'Color','r','LineWidth',3);

%plot(t,y0,'r')

axis tight
line([0,0],ylim,'Color','r','LineWidth',2,'LineStyle','--')
line([20,20],ylim,'Color','r','LineWidth',2,'LineStyle','--')
xlabel('time from stimulus onset [s]')
ylabel('delta cosine distance')
legend(b, labs,'Box','off')
hold off


%%

baseline_tag = cat(2,stims{[1,2]});

secrange = [-5, 40];

data = experiment.series{todo_fish(1)}.data;
fs = data.meta.framerate/data.meta.downsample;
interval = [floor((data.stim_on_sec+secrange(1))*fs) : ...
    1+floor((data.stim_on_sec+secrange(2))*fs)];

allwithintrialcorrintime = nan(length(interval));
withintrialcorrintimeperfish = nan(length(interval));
alldDdt = []

for i_fish = 1:numel(todo_fish)
    
    data = experiment.series{todo_fish(i_fish)}.data;
    fs = data.meta.framerate/data.meta.downsample;
    interval = [floor((data.stim_on_sec+secrange(1))*fs) : ...
        1+floor((data.stim_on_sec+secrange(2))*fs)];
    
    todo_trials = find(~ismember(data.stim_type, baseline_tag));
    for i_trial = 1:numel(todo_trials)
        tmp = data.singletrial{todo_trials(i_trial)}.intime.distancemat{1};
        tmp = tmp(interval,interval);
        allwithintrialcorrintime(:,:,end+1) = ...
            nan(size(allwithintrialcorrintime,1),size(allwithintrialcorrintime,2));
        allwithintrialcorrintime(1:length(interval),1:length(interval),end) = tmp;
        
        tmp = data.singletrial{todo_trials(i_trial)}.intime.ddistancedt;
        tmp = tmp(interval);
        alldDdt = [alldDdt; nan(1,size(alldDdt,2))];
        alldDdt(end,1:length(interval)) = tmp;
    end
    
%     withintrialcorrintimeperfish(:,:,i_fish) = ...
%         nan(size(withintrialcorrintimeperfish,1),size(withintrialcorrintimeperfish,2));
%     withintrialcorrintimeperfish(:,:,i_fish) = ...
%         nanmean(allwithintrialcorrintime(:,:,end-numel(todo_trials)+1:end),3);

end
allwithintrialcorrintime = allwithintrialcorrintime(:,:,2:end);
withintrialcorrintimeperfish = withintrialcorrintimeperfish(:,:,2:end);


t = linspace(secrange(1),secrange(2),size(allwithintrialcorrintime,1));
% y = 1-nanmean(allwithintrialcorrintime,3);
y = 1-allwithintrialcorrintime(:,:,5);
figure; imagesc(t(1:end-2),t(1:end-2),y(1:end-2,1:end-2))
caxis([.3, 1])
axis square
xlabel('time from stimulus onset [s]')
ylabel('time from stimulus onset [s]')


figure; hold on
y = nanmean(alldDdt,1); 
y = y - nanmean(y(1:4));
curve1 = y + std(alldDdt,[],1,'omitnan'); curve1 = curve1(1:end-2);
curve2 = y - std(alldDdt,[],1,'omitnan'); curve2 = curve2(1:end-2);
y = y(1:end-2);
h = patch([t(1:end-2),fliplr(t(1:end-2))],[curve1, fliplr(curve2)],'k','FaceAlpha',.3,'EdgeColor','none')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t(1:end-2), y,'k','LineWidth',3)
xlabel('time from stimulus onset [s]')
ylabel('dXi/dt')
axis tight
range = abs(nanmin(y) - nanmax(y));
% ylim([nanmin(y)-.1*range, nanmax(y)+.1*range])
view([90 -90])
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
box off; hold off

set(gcf, 'Position', [50 50 200 300]);























