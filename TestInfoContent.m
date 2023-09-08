% This function follows a bayesian approach and asks: what is the
% probability for the stimulus identity S to be equal to s (one particular
% identity), based on the average activity vector of the network in a small
% time bin?
% In PDF mode, it first extracts the PDF of any neuron's activity as the
% normalized histogram over all trials and timeframes of the defined "stimulus window"
% (may be different from the actual stimulus window). The independent
% probability P(Xi) of any activity vector is then the product of the independent
% probailities for each neuron's activity value, computed as the value
% associated with the PDF bin in which it is contained.
% Then, the function computes the PDF of each neuron's activity
% analogously, but only considering trials in which that stimulus type was
% presented. The resulting P(Xi|S=s) is the conditional probability of the
% average activity vector in the considered time bin, given the stimulus
% identity of interest.
% Finally, the function computes the conditional probability that S=s,
% given the observation of the average activity vector Xi in the current
% time bin.





traces = data.tracesdns;
L = data.Ldns;
N = data.N;
nbins = 50;
pl = true;
method = 'PDF';


interval_step = 3;
interval_len = 5;
interval = 1:(interval_len*fs+1);

stimintvstart = data.stim_on_sec;
stimintvend = data.stim_off_sec+11;

allstimtypes = unique(data.stim_type);



hf = {};

for i_stim = 1:numel(allstimtypes)
    thisstimtype = allstimtypes(i_stim);
    
    
% get P(S=s)
ntrials = numel(data.trials);
thisstimtrials = ismember(data.stim_type,thisstimtype);
PS = sum(thisstimtrials)/ntrials;


% initialize singletrial segment
fs = data.meta.framerate/data.meta.downsample;
sortorder=1:24;
stiminterval = [floor(stimintvstart*fs) : 1+floor(stimintvend*fs)];
stimivtraces = selectTimeFrame(traces,stiminterval,L);

% define PDF(n)
[V, Edges] = PDFfromHistogram(stimivtraces, length(stiminterval), N, nbins);
if pl; plotPDF(V,Edges); end

[partialV, partialEdges] = ...
    PDFfromHistogram(stimivtraces(:,thisstimtrials), length(stiminterval), N, nbins);
if pl; plotPDF(partialV,partialEdges); end


alltP.v = [];
alltP.exp = [];

for i_t = 1:floor((length(stiminterval)-length(interval))/(interval_step*fs))
    
    thisinterval = stiminterval(1)-1 + interval+floor((i_t-1)*interval_step*fs)


    p.v = []; p.exp = [];

    x435 = []; y435 = [];
    %figure(4524); hold on

    for i_trial = 1:ntrials

        % get Xi_avg for this trial
        disp(['trial #',num2str(data.trial_num(i_trial))])
        win_all = selectTimeFrame(traces(:,i_trial), thisinterval, L);
        Xi_avg = regionnanmean(win_all, length(thisinterval)); % avg actvect

        % get P(Xi_avg) considering all available data
        PXi = getPXi(Xi_avg, V, Edges, stimivtraces, length(stiminterval), N, method);
        [v, exp] = preciseProd(PXi, N);
        PXi_tot.v = v; PXi_tot.exp = exp;
        % now P(Xi) = v*10^exp


        % get P(Xi_avg | S=s)
        PXigivenS = getPXi(Xi_avg, partialV, partialEdges, ...
            stimivtraces(:,thisstimtrials), length(stiminterval), N, method);
        [v, exp] = preciseProd(PXigivenS, N);
        PXigS_tot.v = v; PXigS_tot.exp = exp;



        % get P(S=s | Xi)
        tmpv = PXigS_tot.v * PS / PXi_tot.v;
        tmpexp = PXigS_tot.exp - PXi_tot.exp;
        PSgXi.v = tmpv; PSgXi.exp = tmpexp;
        
%         if tmpexp>0; dbstop; end
        p.v = [p.v; tmpv]; 
        p.exp = [p.exp; tmpexp];



        y = (PXigivenS*PS)./PXi;
        
%         figure(4524)
%         x = .2*randn(1,N) + ones(1,N)*i_trial;
%         if thisstimtrials(i_trial)
%             scatter(x, y, 5, 'g', 'filled')
%         else
%             scatter(x, y, 5, 'b', 'filled')
%         end
%         line(i_trial+[-.3,.3], mean(y)*[1 1], 'Color', 'm', 'LineWidth', 2)
%         line(i_trial+[-.2,.2], median(y)*[1 1], 'Color', 'r', 'LineWidth', 2)


        x435 = [x435, 1+.05*randn(1) + thisstimtrials(i_trial)*1];
        y435 = [y435, mean(y)-median(y)];
    end


%     figure;hold on
%     xlim([0 3]); xticks([1 2]); xticklabels({'other odors', 'thisodor'})
%     scatter(x435,y435,20,'b','filled')
%     boxplot(y435,thisstimtrials)
%     pval1 = kruskalwallis(y435,thisstimtrials,'off');
%     title(['p: ',num2str(pval1)])

        
    %figure;
%     semilogy(p.v.*10.^p.exp)
%     hold on
%     scatter(find(thisstimtrials),p.v(thisstimtrials).*10.^p.exp(thisstimtrials))

%     figure;hold on
%     y = p.exp;
%     scatter(x435,y,20,'b','filled')
%     boxplot(y,thisstimtrials)
%     pval2 = kruskalwallis(y,thisstimtrials,'off');
%     xlim([0.5 2.5]); xticks([1 2]); xticklabels({'other odors', 'thisodor'})
%     title(['joint p: ',num2str(pval2)])

    

    alltP.v = [alltP.v p.v];
    alltP.exp = [alltP.exp p.exp];
    
end

%% common plot

hf{end+1} = figure; hold on
[pval, medsthis, medsnotthis] = deal([]);
title(thisstimtype{1})
xlabel('time bin #')
ylabel('prob. this stimulus')
xlim([0 size(alltP.v,2)+1])
xlabs = {};

for i_t = 1:size(alltP.v,2)
    
    ythis = alltP.v(thisstimtrials,i_t).*10.^alltP.exp(thisstimtrials,i_t);
    ynotthis = alltP.v(~thisstimtrials,i_t).*10.^alltP.exp(~thisstimtrials,i_t);
    
    medsthis = [medsthis; median(ythis)];
    medsnotthis = [medsnotthis; median(ynotthis)];
    
    x = .05*randn(1,numel(ythis)) + ones(1,numel(ythis))*i_t;
    scatter(x, ythis, 10, 'r', 'filled')
    
    x = .05*randn(1,numel(ynotthis)) + ones(1,numel(ynotthis))*i_t;
    scatter(x, ynotthis, 10, 'b', 'filled')
    
    pvaltmp = kruskalwallis(alltP.v(:,i_t).*10.^alltP.exp(:,i_t), ...
        thisstimtrials,'off');
    
    if pvaltmp<.05
        a = scatter(i_t, 100*max([median(ythis), median(ynotthis)]), 150, ...
            'Marker', '*', 'MarkerEdgeColor', 'k');
    end
    
    pval = [pval, pvaltmp];
    
    thisinterval = [stiminterval(1)-1 + interval+floor((i_t-1)*interval_step*fs)]/fs;
    xlabs{end+1} = [num2str(round(thisinterval(1))),'-', ...
        num2str(round(thisinterval(end))),'s'];

end

b(1) = plot(medsthis, 'k-', 'LineWidth', 2);
b(2) = plot(medsnotthis, 'k:', 'LineWidth', 2);
if exist('a','var') && ~isempty(a); b(3) = a; end
legend(b, {'this stim', 'other stim', 'significance'})

xticks(1:size(alltP.v,2)); xticklabels(xlabs); xtickangle(45)

set(gca, 'YScale', 'log')

end


%% save plots

fpathout = [FileIn_path,'\',Series_ID,'\figures\',mfilename];
if ~exist(fpathout,'dir'); mkdir(fpathout);end

for i_f = 1:numel(hf)
    
    nameflag = ['identity_',allstimtypes{i_f}]
    filename = [fpathout,'\',nameflag];
    saveas(hf{i_f},[filename,'.jpeg'])
    savefig(hf{i_f},filename)
    
end




%% functions


function pXi = getPXi(Xi_avg, V, Edges, thistraces, L, N, method)

pXi = ones(N,1);
for i_n = 1:N
    switch method
        case 'PDF'
            tmp = pThisValuefromPDF(Xi_avg(i_n), V(i_n,:), Edges(i_n,:));
        case 'percentile'
            x = L*(i_n-1)+1:L*i_n;
            thisneuron = selectTimeFrame(thistraces,x);
            tmp = calcPercentile(thisneuron(:), Xi_avg(i_n));
    end
    pXi(i_n) = tmp;
end

end


function x = calcPercentile(data, val)
    
nless = sum(data<val);
nequal = sum(data==val);
x = (nless + .5*nequal)/length(data);

end


function [v, exp] = preciseProd(pXi, N)

%alright need to get around double precision... here it goes
    exp = 0; v = 1;
    for i_n = 1:N
        if pXi(i_n)==0; continue;end
        thisexp = floor(log10(pXi(i_n)));
        v = v*(pXi(i_n)/10^(thisexp+1));
        exp = thisexp+1 + exp + floor(log10(v))+1;
        v = v/10^(floor(log10(v))+1);
    end
    
end


function p = pThisValuefromPDF(actval, thisV, thisEdges)

%identify correponding PDF coordinate
[~,idx] = max(thisEdges( thisEdges <= actval )); 
if idx==length(thisEdges) ; idx = idx-1; end

% extract probability to see this value
tmp = thisV(idx);
if tmp==0
    % find previous non-zero
    prev = idx-1;
    while thisV(prev) == 0; prev = prev-1; end
    % find next non-zero
    next = idx+1;
    while thisV(next) == 0; next = next+1; end
    %set probability as interp value
    span = next-prev;
    tmp = (idx-prev) * ((thisV(next)-thisV(prev))/span) + thisV(prev);
end
if isempty(tmp); tmp = 0; end

p = tmp;

end


function [V, Edges] = PDFfromHistogram(traces, L, N, nbins)

V = zeros(N,nbins);
Edges = zeros(N,nbins+1);
for i_n = 1:N
    x = L*(i_n-1)+1:L*i_n;
    thisneuron = selectTimeFrame(traces,x);
    
    [v, e] = histcounts(thisneuron(:), nbins);
    V(i_n,:) = v / numel(thisneuron(:));
    Edges(i_n,:) = e;
end

end


function plotPDF(V,Edges)

figure;
subplot(121); imagesc(V)
xlabel('bin #'); ylabel('neuron #')
title('normalized activity level histogram'); colorbar
subplot(122); imagesc(Edges)
xlabel('bin #'); ylabel('neuron #')
title('histogram edges'); colorbar

end









