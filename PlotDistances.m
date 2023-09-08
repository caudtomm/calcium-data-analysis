

%% distances across trials by trial number

% whole trial
distTraceData(data,[],[],[],'cosine',cells); % cosine distance
distTraceData(data,[],[],[],'correlation',cells); % 1-pearson

% response window only
interval = [floor(data.stim_on_sec*fs) : 1+floor(data.stim_off_sec*fs)];
distTraceData(data,[],[],[],'cosine',cells,interval); % cosine distance
distTraceData(data,[],[],[],'correlation',cells,interval); % 1-pearson

%% response window only + top responsive cells

FindTopUnits
interval = [floor(data.stim_on_sec*fs) : 1+floor(data.stim_off_sec*fs)];
distTraceData(data,[],[],[],'cosine',alltopunits,interval); % cosine distance
distTraceData(data,[],[],[],'correlation',alltopunits,interval); % 1-pearson


%% distances across trials by stim type

idx = data.idx_by_stim_type;

% whole trial
distTraceData(data, idx); % cosine distance
distTraceData(data,idx,[],[],'correlation'); % 1-pearson

% response window only
interval = [floor(data.stim_on_sec*fs) : 1+floor(data.stim_off_sec*fs)];
distTraceData(data,idx,[],[],'cosine',cells,interval); % cosine distance
distTraceData(data,idx,[],[],'correlation',cells,interval); % 1-pearson


%% cosine dist across trials using sliding time window

cax = [.1 1];
method = 'cosine';

[~,idx] = sort(data.trial_num);
[avgd, semd, stdd] = plotCosDMatInTime(data,2,1,cells,idx,1,cax,'trialnum',method);
close all

%% more dist across trials using sliding time window

method = 'cosine';

idx = data.idx_by_stim_type;
plotCosDMatInTime(data,5,3,cells,idx,1,cax,'stimtype',method)
close all

cax = [.1 1];
method = 'correlation';

idx = sort(data.trial_num);
plotCosDMatInTime(data,5,3,cells,idx,1,cax,'trialnum',method)
close all

idx = data.idx_by_stim_type;
plotCosDMatInTime(data,5,3,cells,idx,1,cax,'stimtype',method)
close all

%% distances across time per trial

pathout = 'figures';
if ~isfolder(pathout); mkdir(pathout); end
cd(pathout)

data = singleTrialCorrInTime(data,data.stim_on_sec,data.stim_off_sec,1,1)
fprintf('Saving...'); save(DATAin,'data'); fprintf(' DONE!\n')

CountResponseStages; % outputs figure FIG and JPEG to pwd

cd(FileIn_path)

%% distances across time and across trials

plotFishDistInTime(data,[-5,15],'cosine',1);

disp('DONE!')

