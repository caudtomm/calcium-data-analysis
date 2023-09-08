
%% postprocessing incipit

% load data struct
DATAin = strcat(FileIn_path,'\',Series_ID,'_DATA.mat');
fprintf(strcat('\n\nLoading: ',strrep(DATAin,'\','\\'),' ...'))
load(DATAin); fprintf(' DONE!\n\n')

% show stuff on console
disp(strcat('Period [frames]: ',num2str(data.L)));
disp(strcat('Period [seconds]: ',num2str(data.L/data.meta.framerate)));
disp(strcat('Number of neurons: ',num2str(data.N)));
disp(strcat('Number of trials: ',num2str(numel(data.trials))));
disp(' ')

fs = data.meta.framerate;

data.anat_regions.names


%% select cells in anatomical regions of interest

% regions = data.anat_regions.names;
regions = {'pDp'};

cells = [];
for i_reg = 1:numel(regions)
    [~,idx] = ismember(regions{i_reg},data.anat_regions.names);
    cells = [cells; data.anat_regions.cells{idx}];
end
cells = sort(unique(cells));
[~,cells] = ismember(cells,data.common_units);
cells(cells==0) = [];


