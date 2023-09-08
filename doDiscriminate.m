function doDiscriminate(experiment, classifier, dozscore, stims2use, trainblockmode, s)

global tframe 

nfish = numel(experiment.series);
todo_fish = [1:5];
if isempty(stims2use); stims2use = {'Trp','Ser','Ala','Food'}; end
nstims = numel(stims2use);
trialn2usein = 1:5;
if isempty(s); s = 0; end

oneblock = [1:nstims];
switch trainblockmode
    case '3blocks'
        trainblocksets = {1:3, 2:4, 3:5};
    case 'single'
        trainblocksets = {1,2,3,4,5};
    otherwise
        error('specified training blocks mode is unknown')
end

nsets = numel(trainblocksets);

nshuffles = 50;

if isempty(classifier); classifier = "SVM"; end

% time period to include (relative to [stim_on, stim_on) [s]
if ~exist('tframe','var') || isempty(tframe); tframe = [.5 20]; end


%% action
correctLab = [];
correctTest = [];
correctLabSH = [];
correctTestSH = [];
for i_fish = 1:nfish
    if ~ismember(i_fish,todo_fish);continue;end
    [a, X, stims] = selectData(experiment, i_fish,stims2use,trialn2usein);
%     if dozscore; a = transpose(zscore(a')); end
    
    ncells = size(a,2);
    ntrials = numel(X);
    L = size(a,1)/ntrials;
    
%     a  = zscore(a);
    a = traceFormat(a,L);
    a = squeeze(nanmean(a,1));
    if dozscore; a = zscore(a); end

%     toskipidx = strfind(X,[3 3]); toskipidx = [toskipidx, max(toskipidx+1)]; % idiosyncratic
    
    for i_set = 1:nsets
%         trials_train = reshape([(trainblocksets{i_set}'-1)*max(oneblock)+oneblock]', ...
%             [numel(oneblock)*numel(trainblocksets{i_set}),1]);
        
        trials_train = [];
        for i_stim = 1:nstims
            thisstim = find(ismember(stims,stims2use(i_stim)));
            thisstimidx = find(ismember(X,thisstim));
            if isempty(thisstimidx); continue;end
            trials_train = [trials_train; thisstimidx(trainblocksets{i_set})];
        end

%         trials_train(trials_train >= min(toskipidx)) = ...
%             max(toskipidx) - min(toskipidx) + ...
%             trials_train(trials_train >= min(toskipidx));
        trials_test = [1:ntrials]'; trials_test(trials_train) = [];
    
        trainlabs = stims(X(trials_train))';
        testlabs = stims(X(trials_test))';

        for i = 1:nshuffles+1 % 1x data + 50x shuffle
            switch i
                case 1
                    tmp = a;
                otherwise
                    a_shuf = a;
                    for i_trial = 1:ntrials
                        a_shuf(i_trial,:) = a_shuf(i_trial,randperm(ncells));
                    end
                    tmp = a_shuf;
            end
%             tmp = tmp - nanmean(tmp,1); % men-subtract each feature (neuron)
            trainData = tmp(trials_train,:);
            testData = tmp(trials_test,:);
    
            switch classifier
                case "SVM"
                    [svm, accuracy, predictions] = trainSVM(trainData,trainlabs,stims2use);
                    yfit = svm.predictFcn(testData);
                case "template_match"
                    [templates, templates_weights] = deal(nan(nstims, size(trainData,2)));
                    for i_stim = 1:nstims
                        thisstim_trials = ismember(trainlabs,stims2use{i_stim});
                        templates(i_stim,:) = nanmean(trainData(thisstim_trials,:),1);
%                         tmp = std(trainData(thisstim_trials,:),1);
%                         tmp = nanmin(tmp) ./ tmp;
%                         templates_weights(i_stim,:) = tmp;
                    end

                    method = 'cosine';
                    distances = nan(numel(trials_train),nstims);
                    for i_trial = 1:numel(trials_train)
                        thistrial = trainData(i_trial,:);
                        for i_template = 1:nstims
                            distances(i_trial,i_template) = ...
                                pdist([templates(i_template,:); thistrial],method);
                        end
                    end
                    [distances, idx] = sort(distances,2);
                    predictions = stims2use(idx(:,1))';
                    predictions_confidence = diff(distances(:,1:2),[],2) ./ 2; % DIVISION BY 2 IN THE CASE OF COSINE DISTANCE


                    distances = nan(numel(trials_test),nstims);
                    for i_trial = 1:numel(trials_test)
                        thistrial = testData(i_trial,:);
                        for i_template = 1:nstims
                            distances(i_trial,i_template) = ...
                                pdist([templates(i_template,:); thistrial],method);
                        end
                    end
                    [distances, idx] = sort(distances,2);
                    yfit = stims2use(idx(:,1))';
                    yfit_confidence = diff(distances(:,1:2),[],2);
                    
                case "template_match_adaptive"
                    [templates, templates_weights] = deal(nan(nstims, size(trainData,2)));
                    for i_stim = 1:nstims
                        thisstim_trials = ismember(trainlabs,stims2use{i_stim});
                        templates(i_stim,:) = nanmean(trainData(thisstim_trials,:),1);
%                         tmp = std(trainData(thisstim_trials,:),1);
%                         tmp = nanmin(tmp) ./ tmp;
%                         templates_weights(i_stim,:) = tmp;
                    end

                    method = 'cosine';
                    distances = nan(numel(trials_train),nstims);
                    for i_trial = 1:numel(trials_train)
                        thistrial = trainData(i_trial,:);
                        for i_template = 1:nstims
                            distances(i_trial,i_template) = ...
                                pdist([templates(i_template,:); thistrial],method);
                        end
                    end
                    [distances, idx] = sort(distances,2);
                    predictions = stims2use(idx(:,1))';
                    predictions_confidence = diff(distances(:,1:2),[],2);
                    

                    distances = nan(numel(trials_test),nstims);
                    for i_trial = 1:numel(trials_test)
                        thistrial = testData(i_trial,:);
                        for i_template = 1:nstims
                            distances(i_trial,i_template) = ...
                                pdist([templates(i_template,:); thistrial],method);
                        end
                    end
                    [distances, idx] = sort(distances,2);
                    yfit = stims2use(idx(:,1))';
                    yfit_confidence = diff(distances(:,1:2),[],2) ./ 2; % DIVISION BY 2 IN THE CASE OF COSINE DISTANCE

                otherwise
                    error('unknown classifier')
            end
    
            correct_test = cellfun(@isequal, testlabs, yfit);
            correct_train = cellfun(@isequal, trainlabs, predictions);
            thiscorrect = nan(ntrials,1);
            thiscorrect(trials_test) = correct_test;
            thiscorrect(trials_train) = correct_train;

            switch i
                case 1
                    try
                        correctLab = [correctLab, thiscorrect];
                        correctTest = [correctTest, correct_test];
                    catch
%                         correctLab = [correctLab; nan(15,25)];
%                         correctTest = [correctTest; nan(12,25)];
%                         correctLabSH = [correctLabSH; nan(15,1250)];
%                         correctTestSH = [correctTestSH; nan(12,1250)];
                        
                        correctLab = [correctLab; nan(5,25)];
                        correctTest = [correctTest; nan(4,25)];
                        correctLabSH = [correctLabSH; nan(5,1250)];
                        correctTestSH = [correctTestSH; nan(4,1250)];

                        correctLab = [correctLab, thiscorrect];
                        correctTest = [correctTest, correct_test];
                    end
                otherwise
                    correctLabSH = [correctLabSH, thiscorrect];
                    correctTestSH = [correctTestSH, correct_test];
            end
            
            
        end
    end
end

correctLab = permute(traceFormat(correctLab',nsets),[2,3,1]);
correctTest = permute(traceFormat(correctTest',nsets),[2,3,1]);

correctLabSHtmp = permute(traceFormat(correctLabSH',nsets*nshuffles),[2,3,1]);
correctLabSH = nan(size(correctLab));
for i_set = 1:nsets
    idx = (i_set-1)*nshuffles + 1 : i_set*nshuffles;
    correctLabSH(:,:,i_set) = nanmean(correctLabSHtmp(:,:,idx),3);
end
clear correctLabSHtmp

correctTestSHtmp = permute(traceFormat(correctTestSH',nsets*nshuffles),[2,3,1]);
correctTestSH = nan(size(correctTest));
for i_set = 1:nsets
    idx = (i_set-1)*nshuffles + 1 : i_set*nshuffles;
    correctTestSH(:,:,i_set) = nanmean(correctTestSHtmp(:,:,idx),3);
end
clear correctTestSHtmp

%% plotting
figure;set(gcf,'Color','w')
for i_set = 1:nsets
    subplot(2,nsets,i_set); imagesc(correctLab(:,:,i_set))
    title(['training set ',num2str(i_set)])
    ylabel('fish')
end
for i_set = 1:nsets
    subplot(2,nsets,i_set+nsets); imagesc(correctLabSH(:,:,i_set))
end
xticks(1:numel(X)); xticklabels(stims(X));


C = optimalcolors(nsets+1); C = C(2:end,:);

h = figure; hold on
t = 1:ntrials;
for i_set = 1:nsets
    tmp = correctLab(:,:,i_set);
    ymean = 100* nanmean(tmp,1);
    ystd = 100* squeeze(std(tmp,[],1,'omitnan'));

    curve1 = ymean + ystd; curve2 = ymean - ystd;
    hp = patch([t,fliplr(t)],[curve1'; fliplr(curve2)'], ...
        C(i_set,:),'FaceAlpha',.3,'EdgeColor','none');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(t,ymean,'Color',C(i_set,:),'LineWidth',2,'DisplayName',['training set ',num2str(i_set)]);
end
u = legendUnq();
legend(u,'Box','off','color','none','Location','best','EdgeColor','w','TextColor','w')
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
grid off
axis tight
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 700 200]);
xticks(1:ntrials); xticklabels(stims(X));
ylabel('% hits')

pathout = fullfile('W:\scratch\gfriedri\caudtomm\ev_data\experiments\figures\');
FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset','.fig'));
if s; savefig(h,FileOut); end
FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset','.svg'));
if s; plot2svg(FileOut,h); end




totFractionCorrectLab = squeeze(nanmean(correctTest,2));
totFractionCorrectLabSH = squeeze(nanmean(correctTestSH,2));

lip = .4;
h= figure; hold on
ylim([-.1,1.1]), xlim([1-lip, nsets+lip])
csh = [.6 .6 .6];
chancelv = 1/nstims;
% chancelv = 1/4
line(xlim,repelem(chancelv,2), ...
    'Color','c','LineWidth',3,'LineStyle',':','DisplayName','chance')
for i_set = 1:nsets
    x = i_set+[-1, 1]*lip*2/3;
    
    curve1 = repelem(max(ylim),2); curve2 = repelem(min(ylim),2);
    hp = patch([x,fliplr(x)],[curve1'; fliplr(curve2)'], ...
        C(i_set,:),'FaceAlpha',.3,'EdgeColor','none');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';

    line(x,repelem(nanmean(totFractionCorrectLab(:,i_set)),2), ...
        'Color','r','LineWidth',5)
    line(x,repelem(nanmean(totFractionCorrectLabSH(:,i_set)),2), ...
        'Color',csh,'LineWidth',5,'DisplayName','shuffle')
end


jit = .1*randn(size(totFractionCorrectLab));
x = repmat(1:nsets,[size(totFractionCorrectLab,1),1]) + jit;
scatter(x, totFractionCorrectLab,50,'w','filled')
plot(x',totFractionCorrectLab','w','DisplayName','data')

u = legendUnq();
legend(u,'Box','off','color','none','Location','best','EdgeColor','w','TextColor','w')

xticks(1:nsets)
xlabel('training set')
ylabel('overall performance on test data')
grid off
set(gca, 'color', 'none', 'XColor','w', 'YColor','w', 'ZColor','w');
set(gcf, 'color', 'none'); 
set(gcf, 'Position', [50 50 200 400])

FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset2','.fig'));
if s; savefig(h,FileOut); end
FileOut = fullfile(pathout,strcat(classifier,'classified_by_trainset2','.svg'));
if s; plot2svg(FileOut,h); end

[p,h,stats] = ranksum(totFractionCorrectLab(:,1),totFractionCorrectLab(:,end))
y = totFractionCorrectLab(:,end-2:end);
[p,h,stats] = ranksum(totFractionCorrectLab(:,1),y(:))

