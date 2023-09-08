function timetraces = desurround(plane,period,pl)


 %%% center-surround
    ROI_map = plane{1}.ROI_map;
    fs = plane{1}.meta.framerate;
    fulltraces = plane{1}.timetraces_raw;
    if nargin<2; period = [1:size(fulltraces,1)]'; end
    movie = plane{1}.movie; movie = movie(:,:,period);

    traces = fulltraces(period,:);
    timetraces = nan(size(traces));
    
    allrois = ROI_map;
    allrois = allrois > 0;
    allrois = imdilate(allrois,strel('disk',6));
    
    cells = unique(ROI_map(:));
    cells(cells==0) = [];
    
    for i_cell = 1:numel(cells)
        roi = plane{1}.ROI_map;
        roi = roi==cells(i_cell);
        
        % get raw signal from surround
        y = imdilate(roi,strel('disk',20)) - imdilate(roi,strel('disk',6)) - allrois;
        y = y>0; y = double(y); y(y==0) = nan;
        surround = movie .* repmat(y,[1,1,size(movie,3)]);
        surround = squeeze(nanmean(surround,[1,2]));
        
        x = traces(:,i_cell);
        s = x - surround;
        baselinemask = abs(s-nanmedian(s)) < std(s,'omitnan')/2;
        k = nanmean(x(baselinemask)./surround(baselinemask));
        s1 = x - k*surround;
        s0 = nanmean(s1(baselinemask));
        s2 = (s1-s0)/(abs(s0)+1);

        t = 0:1/fs:length(x)/fs-1/fs;
        
        if i_cell<=5 && pl
            figure
            subplot(411); plot(t,x,t,surround)
            legend('native center signal','surround signal')
            axis tight; ax(1) = gca;
            subplot(412); plot(t,s,t,5*double(baselinemask))
            legend('naive center-surround','baseline periods')
            axis tight; ax(2) = gca;
            subplot(413); plot(t,s1,t,s0)
            subplot(413); plot(t,s1,t,s0*ones(1,length(t)))
            legend('dF','baseline intensity')
            axis tight; ax(3) = gca;
            subplot(414); plot(t,s2)
            xlabel('time [s]')
            legend('dF/F')
            axis tight; ax(4) = gca;
            linkaxes(ax,'x')
        end
        
        timetraces(:,i_cell) = s2;
    end