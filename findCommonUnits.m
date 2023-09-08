function [is_common, is_unit] = findCommonUnits(traces, nan_th, inf_th)
% outputs
% is_common : array of indices corresponding to an existing and 
% non-nan timetrace in all source datasets
% is_unit : array of indices corresponding to an existing and 
% non-nan timetrace in each source dataset
%
% inputs
% traces : is a cell array. each item of the array can be either the name 
% of a defROIs() output MAT file (from which the timetraces will be 
% extracted), or directly the timetraces in matrix form (rows are cells).
% if timetrace data can't be found, file is skipped and the result shows
% units common to all remaining data frames.
% nan_th : is a number >=0 and <1. Indicates the maximum allowed percentage
% of NaNs in any one trace for one neuron, in order for the trace to be
% considered valid. nan_th = 0 by default.
%
% 
% Tommaso Caudullo 08.04.2020
% Tommaso Caudullo 21.01.2020

%setup workspace
switch nargin
    case 0
        error('Not enough input arguments!')
    case 1
        if ~iscell(traces); error('Incorrect data type!'); end
        nan_th = 0;
        inf_th = 0;
    case 2
        if ~iscell(traces); error('Incorrect data type!'); end
        if isempty(nan_th); nan_th = 0; end
        if ~isnumeric(nan_th); error('Incorrect data type!'); end
        if nan_th<0 || nan_th>=1; error('Unacceptable value for nan_th input var!'); end
        inf_th = 0;
    case 3
        if ~iscell(traces); error('Incorrect data type!'); end
        if isempty(nan_th); nan_th = 0; end
        if ~isnumeric(nan_th); error('Incorrect data type!'); end
        if nan_th<0 || nan_th>=1; error('Unacceptable value for nan_th input var!'); end
        if isempty(inf_th); inf_th = 0; end
        if ~isnumeric(inf_th); error('Incorrect data type!'); end
        if inf_th<0 || inf_th>=1; error('Unacceptable value for inf_th input var!'); end
    otherwise
        error('Too many input arguments!')
end


% initialize vars
is_common = [];
is_unit = {};


% get indices of cells containing data in each data frame
for i_t = 1:numel(traces)
    fprintf(strcat('Trace id #',num2str(i_t),' ... '))
    thistrace = traces{i_t};
    
    % extract numeric timetrace
    if isstring(thistrace) || ischar(thistrace)
        if ~exist(thistrace,'file'); warning('data not found : skipped.'); continue;end
        fprintf('loading...')
        tmp = load(thistrace);
        fprintf('done. ')
        
        thistrace = tmp.plane{1}.timetraces';
        clear tmp
    elseif ~isnumeric(thistrace)
        warning('incorrect data type: skipped.')
    end
    
    % timetrace is now numeric
    % rows are cells
    
    % take which indices contain sufficient unit activity data in this trace
    tmp = sum(isnan(thistrace),2) <= nan_th*size(thistrace,2) & ...
        sum(isinf(thistrace),2) <= inf_th*size(thistrace,2); % returns boolean array whether row contains sufficient data
    is_unit{i_t} = find(tmp); clear tmp
    fprintf(strcat(num2str(length(is_unit{i_t})),' cells found\n'))
end


% find cells common to all data frames
max_N = 0;
tmp = []; % external initialization
for i_u = 1:numel(is_unit)
    if isempty(is_unit{i_u}); continue;end % skipped frames
    if i_u==1; tmp = is_unit{i_u}; continue; end % internal initialization
    
    tmp = tmp(ismember(tmp,is_unit{i_u}));
    max_N = max([max_N,max(is_unit{i_u})]);
end
is_common = tmp; % store and return

disp('program terminated.')
disp(' ')
disp(strcat('Common found: ',num2str(length(is_common))))
disp(strcat('Tot found: about ',num2str(max_N)))



end