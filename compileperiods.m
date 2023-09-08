function periods = compileperiods(frames,armlength,maxval)
periods = [];
if isempty(frames); return; end
% make periods
jumps = find(diff(frames)>1);
periods = [frames(1);frames(sort([jumps;jumps+1]));frames(end)];
periods = transpose(reshape(periods,[2,length(periods)/2]));
% add arms
periods(:,1) = max([periods(:,1)-armlength,repelem(1,size(periods,1))'],[],2);
periods(:,2) = min([periods(:,2)+armlength,repelem(maxval,size(periods,1))'],[],2);
% merge overlaps
if size(periods,1)<=1; return; end
tomerge = logical([0; periods(2:end,1)-periods(1:end-1,2)<=0]);
tmp = []; thisstart = 1;
for i = 1:size(periods,1)
    if ~tomerge(i)
        thisstart = periods(i,1);
    end
    if length(tomerge)==i || ~tomerge(i+1)
        tmp = [tmp; thisstart, periods(i,2)];
    end
end
periods = tmp;
return
end