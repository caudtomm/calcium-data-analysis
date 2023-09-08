function traces_out = correctTraces(traces, mode)
% 
% 
% 


switch nargin
    case 0
        error('Not enough input arguments!')
    case 1
        if ~iscell(traces); error('Incorrect data type!'); end
        mode = 'native';
    case 2
        if ~iscell(traces); error('Incorrect data type!'); end
        if ~ismember(mode,{'native','lininterp'}); error('Unknown mode!'); end
    otherwise
        error('Too many input arguments!')
end

traces_out = {};

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
    
    % take which indices contain any nans or infs in this trace (broken units)
    tmp = sum(isnan(thistrace),2)-size(thistrace,2) < 0 && ...
        sum(isinf(thistrace),2)-size(thistrace,2) < 0; % returns boolean array whether row contains sufficient data
    is_broken = find(~tmp); is_ok = find(tmp); clear tmp
    
    % correct broken units
    tmp = zeros(size(thistrace));
    idx = 1:size(thistrace,2);
    tmp(is_ok,:) = thistrace(is_ok,:);
    if isempty(is_broken)
        warning('No broken units! Houray!')
        traces_out(i_t) = tmp;
        continue
    end
    for i_u = 1:numel(is_broken)
        thisunittrace = thistrace(is_broken(i_u),:);
        ok_dp =  %identify ok data points
        
        %% to be continued
    end
    
    
end






end