function outputtraces = debubbled_desurround(plane,tag)
%b: [time x cells x trials] trial traces

trialtraces = plane{1}.timetraces_raw;
outputtraces = trialtraces;
th_movmean = 100
periodlen_th = round(plane{1}.meta.framerate*10) % 10 sec is period length threshold

% -----------------

for i_trial = 1:size(trialtraces,3)
%     disp(['... trial ',num2str(i_trial)])
    b = trialtraces(:,:,i_trial);

    c = diff(b,[],1);
    d = nanmean(c,2);
    figure(625); subplot(221); imagesc(b'); title('input'); colormap('gray')
    subplot(223); plot(d/std(d,'omitnan')); axis tight; ylabel('mean dv/dt (in SD units)')
    
    % build periods
    idx = find(isoutlier(d,'movmean',th_movmean));
    if ~isempty(idx)
        idx = [idx(1);idx(1+[find(diff(idx)>periodlen_th)])];
        if idx(1) == 1; idx = idx(2:end); end
    end
    if ~isempty(idx) && idx(end) == size(b,1); idx = idx(1:end-1); end
    periods = [[1;idx+1] [idx;size(b,1)]]

    
    % ----------------
    
%     jref= b(periods(1,1):periods(1,2),:);
%     baselinek = nanmedian(jref,'all')
    j = [];
%     k = [];
    
    for i_g = 1:size(periods,1)
        g = periods(i_g,:);
%         tmp = b(g(1):g(2),:);
%         normk = nanmean(tmp,'all')
%         scalingk = nanmean(std(jref,[],2,'omitnan')) / nanmean(std(tmp,[],2,'omitnan'))
%         tmp = [scalingk*(tmp-normk)+baselinek];

        %%% center-surround
        tmp = desurround2(plane,[g(1):g(2)]',0);
    
        j = [j;tmp];
%         k = [k;thisk];
    end
    
    c = diff(j,[],1);
    d = nanmean(c,2);
    figure(625)
    subplot(222); imagesc(j'); title('result'); colormap('gray')
    subplot(224); plot(d/std(d,'omitnan')); axis tight

%     saveas(gcf,['W:\scratch\gfriedri\caudtomm\ev_data\TC_220209_TC0004_02_sxpDp\desurrounding_results\trial',...
%         tag,'.png'])

    outputtraces(:,:,i_trial) = j;

end