function data = checkTracesMan(data)
% prompts to check traces manually
% data = checkTracesMan(data)


to_del = [];
prompt = 'Next bad unit: [''k'' to stop]';

for i_trace = 1:size(data.traces,2)
    thistraces = reshape(data.traces(:,i_trace), data.L, data.N);
    figure; imagesc(thistraces); title(num2str(i_trace)); xlabel('units')
    
    
    while true
        x = input(prompt,'s');
        y = str2num(x);
        if isempty(y)
            if strcmp(x,'k'); break
            else disp('Input invalid.'); continue;
            end
        end
        
        to_del = [to_del; y];
    end
    
    close(gcf)
end

to_del = unique(to_del);
data.common_units(to_del) = [];

data.traces = takeoutunits(data.traces,to_del,data.N,data.L);


data.N = data.N - numel(to_del);

end



function traces = takeoutunits(traces, to_del, N, L)

mask = ones(N,1); mask(to_del) = 0;
mask = repelem(mask,L);

traces = traces(logical(mask),:);

end