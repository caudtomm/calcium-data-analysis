function [idx_out, sortedtraces, sortedtrials, sortednums, sortedtypes, sortedlocalcorr] ...
    = sortbyStimType(data)
% [idx_4, sortedtraces, sortedtrials, sortednums, sortedtypes, sortedlocalcorr] ...
%     = sortbyStimType(data)

a = unique(data.stim_type);

idx_4 = [];
for i_stim = 1:numel(a)
    idx = find(strcmp(data.stim_type,a{i_stim}));
    [b, idx_2] = sort(data.trial_num(idx));
    idx_3 = idx(idx_2);
    
    c = corrcoef(data.traces(:,idx_3));
    c_1_y(i_stim,1:size(c,1)) = c(1,:);
    idx_4 = [idx_4, idx_3];
end

% output
sortedtraces = data.traces(:,idx_4);
sortedtrials = data.trials(idx_4);
sortednums = data.trial_num(idx_4);
sortedtypes = data.stim_type(idx_4);
sortedlocalcorr = data.localCorrelations(:,:,idx_4);

idx_out = idx_4;