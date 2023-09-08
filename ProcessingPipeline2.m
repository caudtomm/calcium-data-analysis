

% tungstenlocation = '\\tungsten-nas.fmi.ch\tungsten';
tungstenlocation = 'W:';
disp('this is a new branch')

% FileIn_path = 'C:\Users\caudtomm\Desktop\for_TC';
% FileIn_path = "W:\scratch\gfriedri\caudtomm\ev_data\TC_210505_TC0004_01_sxpDp";
Series_ID = 'TC_230506_TC0003_230502beh2b2_sxpDp_odorexp004_RPB3144501500AG';
FileIn_path = fullfile(tungstenlocation,"scratch\gfriedri\caudtomm\ev_data",Series_ID);
FileIn_ext = 'tif';
ref = 1;
anatomy = 'anatomy';


auto = 1;
do_prune = 1;
do_register = 0;
save_reg_stacks = 1; % only relevant if do_register

rawtrialsfolder = 'trials';
if do_register
    trialsfolder = rawtrialsfolder;
    clahefolder = 'trials_clahe';
else
%     trialsfolder = 'trials_warp';
    trialsfolder = 'reg_stacks2_raw';
    clahefolder = 'reg_stacks2_clahe';
end

% startup
cd(FileIn_path)

% odor delay to nostril
odordelay = 0; % sec



%% initialize preprocessing

FileIn = strcat(FileIn_path,'\',Series_ID,'_check.mat');
if exist('to_keep','var'); warning('to_keep variable found in workspace: deleted'); end
if exist('to_discard','var'); warning('to_discard variable found in workspace: deleted'); end
if exist('is_complete','var'); warning('is_complete variable found in workspace: deleted'); end
clear to_keep to_discard is_complete
if exist(FileIn,'file')
    load(FileIn);
    if ~exist('to_keep','var'); warning('to_keep variable not found and could not be imported!'); end
    if ~exist('to_discard','var'); warning('to_discard variable not found and could not be imported!'); end
    if ~exist('is_complete','var'); warning('is_complete variable not found and could not be imported!'); end
else; [to_keep,to_discard,is_complete] = deal({});
end


% get list of files with the first reference file at the top
fprintf(strcat('\nFolder: ',strrep(FileIn_path,'\','\\'),'.\n\n'))
files = dir(fullfile(FileIn_path,trialsfolder,strcat(Series_ID,'*.',FileIn_ext)));
n = [];
for i_f = 1:numel(files)
    file = files(i_f);
    fspecs = getFileNameSpecs(file.name);
    n = [n;fspecs.trial_num];
end
n0 = min(n)-1; n=n-n0;
[~,idx] = sort(n); n = n(idx);
files = files(idx);
% files_ref = dir(strcat(FileIn_path,'\',Series_ID,'*',ref,'*.',FileIn_ext));
% if ~isempty(files_ref); files_ref = files_ref(1); end
% file_ref_name = files_ref(1).name;
% tmp = files; n = 0;
% for i = 1:numel(files)
%     if ~isempty(strfind(files(i).name,file_ref_name))
%         tmp(i-n) = []; n = n+1;
%     elseif ~isempty(strfind(files(i).name,anatomy)) % discard anatomy stacks from file list
%         tmp(i-n) = []; n = n+1;
%     end
% end
% files = [files_ref;tmp];
ref_n = find(n==ref);
files_ref = files(ref_n);
files(ref_n) = []; files = [files_ref;files];
clear tmp files_ref

% initialize vars
ref_img = [];
ref_filename = [];


FileIn = strcat(rawtrialsfolder,'\',Series_ID,'_',num2str(n0+1, '%05.f'));
[stim_series, ch2stimtype, ch2stimtype_map] = read_stim_series(FileIn);


if isempty(extractBefore(files(1).name,'_reg'))
    FileIn = fullfile(rawtrialsfolder,files(1).name);
else
    FileIn = fullfile(rawtrialsfolder,strcat(extractBefore(files(1).name,'_reg'),'.tif'));
end
[~,~,fs,~,~,~,~] = read_metadata_function(FileIn);
stim_on = stim_series(1,3) / fs + odordelay;
stim_off = stim_series(1,4) / fs + odordelay;

f0_window = [floor((stim_on-10)*fs) floor(stim_on*fs)-10];
response_window = floor(stim_series(1,3:4) + odordelay*fs);

%% preprocess data

start_f = 1; ref_f = 1;
if numel(is_complete)>=numel(files)
    prompt = 'Trial num to start with: ';
    x = input(prompt, 's');
    x = str2num(x);
    tmp = 1;
    for i_f = 1:numel(files)
        file = files(i_f);
        filename = strcat(file.folder,'\',file.name);
        fspecs = getFileNameSpecs(filename);
        fspecs.trial_num = fspecs.trial_num - n0;
        if fspecs.trial_num == x; tmp = i_f; break; end
    end
    start_f = tmp;
    
    
    prompt = 'Trial num to use as reference: ';
    x = input(prompt, 's');
    x = str2num(x);
    tmp = 1;
    for i_f = 1:numel(files)
        file = files(i_f);
        filename = strcat(file.folder,'\',file.name);
        fspecs = getFileNameSpecs(filename);
        fspecs.trial_num = fspecs.trial_num - n0;
        if fspecs.trial_num == x; tmp = i_f; break; end
    end
    ref_f = tmp;
end




FileIn = strcat(rawtrialsfolder,'\',Series_ID,'_badperiods.csv');
clear allbadperiods
allbadperiods.data = [];
if exist(FileIn,'file'); allbadperiods = importdata(FileIn); end
if isstruct(allbadperiods)
    allbadperiods = allbadperiods.data;
else
    allbadperiods = [];
end


% run through file list and define ROIs
for i_f = start_f:numel(files)
    file = files(i_f);
    
    if isempty(extractBefore(file.name,'_reg'))
        filenameMETADATA = fullfile(rawtrialsfolder,file.name);
    else
        filenameMETADATA = fullfile(rawtrialsfolder,strcat(extractBefore(file.name,'_reg'),'.tif'));
    end
    filename = fullfile(trialsfolder,file.name);
    filenameCLAHE = fullfile(clahefolder,file.name);
    if i_f==1
        fprintf(strcat('\n\nAnalizing reference: ''',file.name,'''\n')) 
    else
        fprintf(strcat('\n\nAnalizing: ''',file.name,'''\n')) 
    end
    
    % get specs from file name
    fspecs = getFileNameSpecs(filename);
    
    
    % skip file if already complete and first round not yet complete
    if numel(is_complete)<numel(files)
        if any(strcmp(is_complete,file.name))
            warning(strcat('Trial: ',file.name,' is complete: skipped.'))
            if i_f == 1
                disp('Reference detected: loading reference image...')
                FileIn = strcat(FileIn_path,'\defROIs\',fspecs.fname,'_defROIs.mat');
                ref = load(FileIn);
                ref_img = ref.plane{1}.anatomy; 
                ref_filename = FileIn;
            end
            continue
        end
    end
    
    
    FileIn = strcat(FileIn_path,'\','defROIs','\',fspecs.fname,'_defROIs.mat');
    part_annotation = [];
    if exist(FileIn,'file') && numel(is_complete)<numel(files)
        part_annotation = FileIn;
    elseif i_f > start_f
%         FileIn = ref_filename;
        FileIn = fullfile('defROIs',strcat(files(i_f-1).name(1:end-length(FileIn_ext)-1),'_defROIs.mat'));
        if exist(FileIn,'file')
            part_annotation = FileIn;
        else
            warning('partial annotation from previous file not found.')
        end
    elseif i_f==start_f && numel(is_complete)>=numel(files)
        % if the first round has been completed, always take the previous
        % item's annotation as partial annotation
        FileIn = fullfile('defROIs',strcat(files(ref_f).name(1:end-length(FileIn_ext)-1),'_defROIs.mat'));
        if exist(FileIn,'file')
            part_annotation = FileIn;
        else
            warning('partial annotation from last file not found.')
        end
    end
    
    shifts = [];
    FileIn = strcat(FileIn_path,'\',fspecs.fname,'_customalign.mat');
    if exist(FileIn,'file')
        disp('Using custom alignment from file...')
        load(FileIn)
    end
    
    fspecs = getFileNameSpecs(filename);
    fspecs.trial_num = fspecs.trial_num - n0;
    fspecs.stim_type = ch2stimtype{fspecs.trial_num};
    
    badperiods = [];
    if ~isempty(allbadperiods); badperiods = allbadperiods(allbadperiods(:,1)-n0==fspecs.trial_num,:); end

    ref_img = strcat(FileIn_path,'\trials_clahe_warpref\',Series_ID,'_commonslice.tif');
    ref_img = double(loadTiffStack(char(ref_img)));
    
    % ----------------------------------------------------- % 
    plane = defROIs(filenameMETADATA,filename,filenameCLAHE,part_annotation,0,do_prune,ref_img,f0_window,response_window,save_reg_stacks,auto,shifts,badperiods,fspecs,do_register);
    % ----------------------------------------------------- % 
    
    %%% time region selective dF/F computation, with center-surround
    %%% calculation
%     if ~do_register; plane{1}.timetraces = debubbled_desurround(plane); end
    


    %%%
    if i_f == start_f
        disp('Reference detected: loading reference image...')
        FileIn = strcat(FileIn_path,'\',fspecs.fname,'_defROIs.mat');
        ref_img = plane{1}.anatomy;
        ref_filename = FileIn;
    end
    
    % file specs are not updated, only stored if not already present in
    % metadata
    FileOut = strcat(FileIn_path,'\defROIs\',fspecs.fname,'_defROIs.mat');
    if ~isfield(plane{1}.meta,'fspecs')
        plane{1}.meta.fspecs = fspecs;
        if exist(FileOut,'file')
            if auto; overwrite = true; else; overwrite = askBoolean('A previous version was found. Overwrite? [y,n] '); end
            if overwrite; save(FileOut,'plane','-v7.3'); else; warning('Changes discarded.'); end
        else
            save(FileOut,'plane','-v7.3');
        end
    end
    
    % is the annotation complete?
    if auto; thisiscomplete = true; else; thisiscomplete = askBoolean('Annotation complete? [y,n] '); end
    if thisiscomplete && ~any(strcmp(is_complete,file.name)); is_complete{end+1} = file.name; end
    
    % keep trial?
    FileIn = strcat(fspecs.fname,'_defROIs.mat');
    if auto; isgood = true; else; isgood = askBoolean('Trial good? [y,n] '); end
    if isgood && ~any(strcmp(to_keep,FileIn))
        to_keep{end+1} = FileIn;
    elseif ~isgood && ~any(strcmp(to_discard,FileIn))
        to_discard{end+1} = FileIn;
    end

    
    % close figures
    if ishandle(88); close(88); end
    if ishandle(199); close(199); end
    
    % keep going to the next trial?
    if auto; go_on = true; else; go_on = askBoolean('Wanna continue to the next trial? [y,n] '); end
    
    % in case you wanted to finish
    if ~go_on; break; end
end
    
% save general info on group
to_keep = unique(to_keep,'stable');
to_discard = unique(to_discard,'stable');
is_complete = unique(is_complete,'stable');
FileOut = strcat(FileIn_path,'\',Series_ID,'_check.mat');
save(FileOut, 'to_keep', 'to_discard', 'is_complete')

%% Update timetraces through debubble + desurround

% k = [];
for i_f = 1:numel(files)
    file = files(i_f);
    file = file.name;
    FileIn = strcat('defROIs\',file(1:end-4),'_defROIs.mat');
    if ~exist(FileIn, 'file'); continue; end
    disp(FileIn)
    load(FileIn);

    % debubbled desurround
    v = debubbled_desurround(plane,num2str(i_f));
    plane{1}.timetraces = v;
%     k = [k;thisk];

    save(FileIn,'plane','-v7.3');

end




%% save group data (common cells only)

% retrieve check file
FileIn = strcat(FileIn_path,'\',Series_ID,'_check.mat');
if exist('to_keep','var'); warning('to_keep variable found in workspace: deleted'); end
if exist('to_discard','var'); warning('to_discard variable found in workspace: deleted'); end
if exist('is_complete','var'); warning('is_complete variable found in workspace: deleted'); end
if exist('unique_ROIs','var'); warning('unique_ROIs variable found in workspace: deleted'); end
clear to_keep to_discard is_complete unique_ROIs
if exist(FileIn,'file')
    load(FileIn);
    if ~exist('to_keep','var'); warning('to_keep variable not found and could not be imported!'); end
    if ~exist('unique_ROIs','var')
        warning('unique_ROIs variable not found and could not be imported! Will be created now as a list of all files in the series')
        
        fprintf('Creating variable and saving ... ')
        
        unique_ROIs = to_keep;
        
        % save general info on group
        FileOut = FileIn;
        save(FileOut, 'to_keep', 'to_discard', 'is_complete','unique_ROIs')
        
        fprintf('DONE!\n')
    end
else; error('check file not found for this data group.')
end

disp(' ');disp(' ')

% identify common units
disp('Identifying common units.')
[is_common, ~] = findCommonUnits(strcat('defROIs\',unique_ROIs),.9,.9);
disp(' ')

% load sample data for initialization
disp('Loading sample data for initialization from:')
FileIn = fullfile(FileIn_path,'defROIs',to_keep{1});
disp(FileIn); load(FileIn)

% initialize vars
disp(' ');
data.L = size(plane{1}.timetraces,1);
data.N = numel(is_common);
data.traces = zeros(data.L*data.N,numel(to_keep));
data.common_units = is_common;
%data.trace_correction = trace_correction;
data.trials = to_keep;
[data.trial_num] = deal(zeros(1,numel(to_keep)));
[data.stim_type] = deal({});
data.meta = plane{1}.meta;
data.localCorrelations = zeros(size(plane{1}.localCorrelations,1),size(plane{1}.localCorrelations,2),numel(to_keep));

%ROI_map of only common units
ROI_map = plane{1}.ROI_map;
isnt_common = setdiff(unique(ROI_map), is_common);
idx = ismember(ROI_map,isnt_common);
ROI_map(idx) = 0;
data.ROI_map_common = ROI_map;

disp(strcat('Period [frames]: ',num2str(data.L)));
disp(strcat('Period [seconds]: ',num2str(floor(data.L/data.meta.framerate))));
disp(strcat('Number of common neurons: ',num2str(data.N)));
disp(strcat('Number of trials: ',num2str(numel(data.trials))));

disp(' ');disp(' ')

for i_f = 1:numel(to_keep)
%     filename = [FileIn_path,'\',to_keep{i_f}];
    filename = fullfile('defROIs',to_keep{i_f});
    disp(strcat('Loading: ',strrep(filename,'\','\\'),' ...'))
    if ~exist(filename,'file'); warning(strcat(strrep(filename,'\','\\'),' not found! - skipped.')); continue; end
    load(filename)
    if ~exist('plane','var'); error('''plane'' variable not found!'); end
    
    common_traces = plane{1}.timetraces(:,is_common(ismember(is_common,[1:size(plane{1}.timetraces,2)])));
    if size(common_traces,1) > data.L
        warning(strcat('Period longer than expected: cut to ',num2str(data.L),' frames.'))
        common_traces = common_traces(1:data.L,:);
    elseif size(common_traces,1) < data.L
        warning(strcat('Period of ',num2str(size(common_traces,1)),' frames is too short! file skipped.'))
    end
    data.traces(:,i_f) = common_traces(:);
    try
        data.trial_num(i_f) = plane{1}.meta.fspecs.trial_num;
    catch ME
        data.trial_num(i_f) = str2double(plane{1}.meta.fspecs.trial_num);
    end
    data.stim_type{i_f} = plane{1}.meta.fspecs.stim_type;
    data.localCorrelations(:,:,i_f) = plane{1}.localCorrelations;
end


% sort data by trial number
[data.trial_num, idx] = sort(data.trial_num);
data.traces = data.traces(:,idx);
data.stim_type = data.stim_type(idx);
data.trials = data.trials(idx);
data.localCorrelations = data.localCorrelations(:,:,idx);
[data.idx_by_stim_type, ~] = sortbyStimType(data);

% trace quality check
data = checkTracesMan(data);

%  denoise cell traces
FileIn = strcat(Series_ID,'_badperiods.csv'); clear allbadperiods
allbadperiods.data = [];
if exist(FileIn,'file'); allbadperiods = importdata(FileIn); end
allbadperiods = allbadperiods.data;
[data.tracesdn, fcutoff] = denoiseCellTraceData(data,0,allbadperiods); % butterworth filt 4 poles

stim_on = stim_series(1,3) / fs + odordelay;
stim_off = stim_series(1,4) / fs + odordelay;

f0_window = [floor((stim_on-10)*fs) floor(stim_on*fs)-10];

% downsample denoised traces
data.meta.downsample = floor(fcutoff*data.meta.framerate/2);
data.tracesdns = downsample(traceFormat(data.tracesdn,data.L),data.meta.downsample);
data.Ldns = size(data.tracesdns,1);
data.tracesdns = traceFormat(data.tracesdns);
data.meta.seriesid = Series_ID;
data.stim_on_sec = stim_on;
data.stim_off_sec = stim_off;
data.f0_window = f0_window;

% select anatomical regions
data = selectAnatRegions(data,false);

disp('... DONE!'); disp(' ')
data

FileOut = strcat(FileIn_path,'\',Series_ID,'_DATA.mat');
fprintf(strcat('Saving: ',strrep(FileOut,'\','\\'),' ...'))
save(FileOut,'data')
fprintf(' DONE!\n')


