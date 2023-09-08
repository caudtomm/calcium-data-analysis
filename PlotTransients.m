%% transients stuff
% OUTDATED

load([FileIn_path '\TC_200108_TC0004_01_sxpDp_noodor_RPB301AG_004_defROIs.mat'])
[~, ~, ~, ~, baseline] = getTransients(plane{1}.timetraces);

%
stim_on = 80;
stim_off = 240;

%
FileIn = [FileIn_path '\' Series_ID '_DATA.mat'];
load(FileIn)
is_common = data.common_units;


% show stuff on console
disp(strcat('Period [frames]: ',num2str(data.L)));
disp(strcat('Period [seconds]: ',num2str(floor(data.L/data.meta.framerate))));
disp(strcat('Number of neurons: ',num2str(data.N)));
disp(strcat('Number of trials: ',num2str(numel(data.trials))));
disp(' ')

% sort data by trial number
[data.trial_num, idx] = sort(data.trial_num);
data.traces = data.traces(:,idx);
data.stim_type = data.stim_type(idx);
data.trials = data.trials(idx);

%
ensembles = [];
trial_nums = [];

%
%files = dir(strcat(FileIn_path,'\',Series_ID,'*',FileIn_ext));
disp(' ')
all_responsive_n_ids = [];
for i = 1:numel(data.trials)
    tmp = getFileNameSpecs(data.trials{i});
    fname = [tmp.fname,'.',tmp.orig_fext];
    disp(fname)
    
    %
    filename = strcat(FileIn_path,'\',fname);
    load(filename)
    
    %
    [activity, on_response_ids, off_response_ids, responsive_n_ids, baseline] = ...
        getTransients(plane{1}.timetraces(:,is_common),baseline,stim_on,stim_off,[],0,0);
    
    %
    tmp = zeros(size(is_common)); tmp(responsive_n_ids) = 1;
    ensembles = [ensembles tmp];
    trial_nums = [trial_nums; plane{1}.meta.fspecs.trial_num];
    all_responsive_n_ids = [all_responsive_n_ids; responsive_n_ids];
    
end
all_responsive_n_ids = unique(all_responsive_n_ids);
disp(strcat('Tot responsive neurons: ',num2str(numel(all_responsive_n_ids)),' of ',num2str(data.N)))

a = unique(data.stim_type);
% x= reshape(trial_nums,[3,6]);
% [~,idx] = sort(x(1,:));
% idx = (idx-1)*3+1;
% idx= repmat(idx,3,1);
% idx(2,:) = idx(2,:)+1;
% idx(3,:) = idx(3,:)+2;

%only responsive neurons
idx = 1:size(ensembles,2);
c = corrcoef(ensembles(all_responsive_n_ids,idx)); % sorted by trial num
figure;imagesc(c); hold on
xticks(1:size(c,1)); xticklabels(trial_nums(idx(:))); xtickangle(90)
yticks(1:size(c,1)); yticklabels(data.stim_type(idx(:)))
title('Correlation of responsive ensembles')
caxis([0 1]); colorbar
hold off

% num of responsive neurons
idx = 1:size(ensembles,2);
c = corrcoef(ensembles(:,idx)); % sorted by stim type
x = sum(ensembles(:,idx),1); % sorted by stim type
figure; bar(x,'m'); hold on
% xticks(1:size(c,1)); xticklabels(repelem(a,3)); xtickangle(90)
xticks(1:size(c,1)); xticklabels(trial_nums(idx(:))); xtickangle(90)
ylabel('n of responsive neurons')


%only responsive neurons
idx = 1:size(ensembles,2);
c = corrcoef(ensembles(all_responsive_n_ids,idx)); % sorted by stim type
c1 = c(1:size(c,1)-1,2:size(c,1));
c2 = c1([(0:size(c,1)-2)*size(c,1)+1]);
figure; bar([0 c2],'FaceColor','c'); hold on;
lab = strcat(num2str(trial_nums(idx(2:end))), ...
    '-',num2str(trial_nums(idx(1:end-1))));
xticks(1:size(c,1)); xticklabels(cat(1,'     ',lab)); xtickangle(90)
ylabel('correlation btw. trials')


%only responsive neurons
% [~, idx] = sort(trial_nums(:));
c = corrcoef(ensembles(all_responsive_n_ids,idx_4)); % sorted by stim type
figure; imagesc(c); hold on
xticks(1:size(c,1)); xticklabels(trial_nums(idx_4(:))); xtickangle(90)
yticks(1:size(c,1)); yticklabels(data.stim_type(idx_4))
title('Correlation of responsive ensembles')
caxis([0 1]); colorbar
hold off


% num of responsive neurons
x = sum(ensembles(:,idx_4),1); % sorted by trial num
figure; bar(x,'m'); hold on
xticks(1:size(c,1)); xticklabels(trial_nums(idx_4(:))); xtickangle(90)
ylabel('n of responsive neurons')


