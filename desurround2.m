function [timetraces] = desurround2(plane,period,pl)


offset = 10;

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

[S_center_M, S_surround_M] = deal(nan(length(period),numel(cells)));
baselinemask_surround_M = false(length(period),numel(cells));
baselinemask_center_M = false(length(period),numel(cells));


for i_cell = 1:numel(cells)
    roi = plane{1}.ROI_map;
    roi = roi==cells(i_cell);
    
    % get raw signal from surround
    y = imdilate(roi,strel('disk',20)) - imdilate(roi,strel('disk',6)) - allrois;
    y = y>0; y = double(y); y(y==0) = nan;
    surround = movie .* repmat(y,[1,1,size(movie,3)]);
    surround = squeeze(nanmean(surround,[1,2]));
    
    
    x = traces(:,i_cell);

    % get baseline periods
%     s = x - surround; baselinemask_center = abs(s-nanmedian(s)) < std(s,'omitnan')/2;
%     s = surround - x; baselinemask_surround = abs(s-nanmedian(s)) < std(s,'omitnan')/2;
%     baselinemask_intersection = baselinemask_surround & baselinemask_center;

    % calc F0
%     F0_center = nanmean(x(baselinemask_center));
%     F0_surround = nanmean(surround(baselinemask_surround));
    F0_center = quantile(x,.1);
    F0_surround = quantile(surround,.1);

    % calc S(t) ["dirty dF"]
    S_center = x - F0_center;
    S_surround = surround - F0_surround;
    
    % calc cell specific k
    [k,~]=demingRegression(S_center_M,S_surround_M,0,1);
    k = min([k,1]);

    dF_overF0 = (S_center-k*S_surround)/(F0_center+offset);
    
    %store vars to mats
%     S_center_M(:,i_cell) = S_center;
%     S_surround_M(:,i_cell) = S_surround;
%     F0_center_M(i_cell) = F0_center;
%     baselinemask_surround_M(:,i_cell) = baselinemask_surround;
%     baselinemask_center_M(:,i_cell) = baselinemask_center;
    
    timetraces(:,i_cell) = dF_overF0;
end


%% determine relationship between S(t)_center and S(t)_surround

% t = 0:1/fs:length(period)/fs-1/fs;

% surroundhigher = S_surround_M > S_center_M;

% x1 = S_center_M(~baselinemask_surround_M & surroundhigher);
% x2 = S_surround_M(~baselinemask_surround_M & surroundhigher);
% b = polyfit(S_center_M,S_surround_M,1);
% 
% %demingregression
% [k,~]=demingRegression(S_center_M,S_surround_M,0,1);
% k = min([k,1]);

if pl


figure; clear ax
subplot(241); imagesc(S_center_M); title('F_c-F0_c'); c = caxis; ax(1) = gca;
subplot(242); imagesc(S_surround_M); title('F_s-F0_s'); caxis(c); ax(2) = gca;
subplot(243); imagesc(S_center_M-k*S_surround_M); title('dF_c = F_c-F0_ck*(F_s-F0_s)'); c = caxis; ax(3) = gca;
subplot(244); imagesc((S_center_M-k*S_surround_M)./F0_center_M); title('dF_c / F0_c'); c = caxis; ax(4) = gca;
subplot(245); histogram(S_center_M); ylabel('histogram')
subplot(246); histogram(S_surround_M);
subplot(247); histogram(S_center_M-k*S_surround_M);
subplot(248); histogram((S_center_M-k*S_surround_M)./F0_center_M);
linkaxes(ax,'xy')

figure; clear ax
subplot(241); imagesc(S_center_M); title('F_c-F0_c'); c = caxis; ax(1) = gca;
subplot(242); imagesc(S_surround_M); title('F_s-F0_s'); caxis(c); ax(2) = gca;
subplot(243); imagesc(timetraces); title('dF_c = F_c-F0_ck*(F_s-F0_s)'); c = caxis; ax(3) = gca;
subplot(245); histogram(S_center_M); ylabel('histogram')
subplot(246); histogram(S_surround_M);
subplot(247); histogram(timetraces);
linkaxes(ax,'xy')

x = S_surround_M;
y = S_center_M;
figure; scatter(x,y,2,'k','filled')
hold on;
line([-20 35],[-20 35],'Color','r','LineWidth',1.5)
plot(x,k*x,'g','LineWidth',1.5);
xlabel('F_s - F0_s')
ylabel('F_c - F0_c')

figure; scatter(S_center_M(~baselinemask_surround_M),...
    S_surround_M(~baselinemask_surround_M),2,'k','filled')
hold on;
scatter(S_center_M(baselinemask_surround_M),...
    S_surround_M(baselinemask_surround_M),2,'b','filled')
line([-20 35],[-20 35],'Color','r','LineWidth',3)
figure; scatter(S_center_M(baselinemask_surround_M),...
    S_surround_M(baselinemask_surround_M),2,'b','filled')
hold on; scatter(S_center_M(~baselinemask_surround_M),...
    S_surround_M(~baselinemask_surround_M),2,'k','filled')
line([-20 35],[-20 35],'Color','r','LineWidth',1.5)
plot(x1,b(2)+b(1)*x1,'g','LineWidth',1.5)
% plot(x1,b(3)+b(1)*x1.*x1+b(2)*x1,'g','LineWidth',1.5)
xlabel('S\_center')
ylabel('S\_surround')
legend({'during surround baseline','during surround non-baseline','identity','surround is higher'})

figure; h = histogram(S_center_M(:));
hedg = h.BinEdges;
hold on; histogram(S_surround_M(:), hedg)
legend({'S\_center','S\_surround'})

end

