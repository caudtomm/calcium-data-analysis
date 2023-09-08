%% do PC1-3 traj graphics

close all

hf = {};


stim_on = data.stim_on_sec;
stim_off = data.stim_off_sec;

ntrials = numel(data.trials);
fs = data.meta.framerate/data.meta.downsample;
nstimtypes = numel(unique(data.stim_type));
stimtypes = unique(data.stim_type);

c = optimalcolors(nstimtypes);

[~,b] = ismember(stimtypes, data.stim_type);

intv_start = stim_on-5;
intv_end = stim_off+20;
interval0 = floor(intv_start*fs):floor((intv_end)*fs);

hf{end+1} = figure(2865);
hf{end+1} = figure(2866);

a1 = {};
a2 = {};

for i_trial = 1:ntrials
   
    y = data.singletrial{i_trial}.intime.pca;
    if isempty(y); continue; end
    y = y(interval0,:);
    
    %y = movmean(y,3);
    interval = 1:size(y,1);
    
     c1 = repmat([linspace(0,.9,length(interval))]',1,3);
     
    [~, stypeid] = ismember(data.stim_type(i_trial), stimtypes);
    
    for i_fig = 1:2
        
        figure(2865+i_fig-1)
        
        if i_fig==1
            a1{end+1}= plot3(y(interval,i_fig),y(interval,1+i_fig),y(interval,2+i_fig), ...
                    'LineWidth', 2, ...
                    'Color', c(stypeid,:));
        elseif i_fig==2
            a2{end+1}= plot3(y(interval,i_fig),y(interval,1+i_fig),y(interval,2+i_fig), ...
                    'LineWidth', 2, ...
                    'Color', c(stypeid,:));
        end
        
        hold on
        
        scatter3(y(interval,i_fig),y(interval,1+i_fig),y(interval,2+i_fig), ...
                15,c1,'filled')
       
        xlabel(['PC',num2str(i_fig)])
        ylabel(['PC',num2str(1+i_fig)])
        zlabel(['PC',num2str(2+i_fig)])
    end
    
    
end

figure(2865); a1a=[];
for i = 1:nstimtypes; a1a = [a1a;a1{b(i)}]; end
legend(a1a,stimtypes, 'Location','NorthEast')


figure(2865); a2a=[];
for i = 1:nstimtypes; a2a = [a2a;a2{b(i)}]; end
legend(a2a,stimtypes, 'Location','NorthEast')








for i_stim = 1:nstimtypes
    thistrials = find(ismember(data.stim_type, stimtypes(i_stim)));
    
    a1a={};
    for i = 1:numel(thistrials)
        a1a{end+1} = ['trial ',num2str(data.trial_num(thistrials(i)))];
    end
    
    a1 = [];
    a2 = [];
    
    for i_trial = 1:numel(thistrials)
        y = data.singletrial{thistrials(i_trial)}.intime.pca;
        if isempty(y); continue; end
        y = y(interval0,:);

        %y = movmean(y,3);
        interval = 1:size(y,1);
        
        c1 = repmat([linspace(0,.9,length(interval))]',1,3);
        
        
        for i_fig = 1:2
            hf{end+1} = figure(3646*i_fig+i_stim);
        
            if i_fig==1
                a1 = [a1; plot3(y(interval,i_fig),y(interval,1+i_fig),y(interval,2+i_fig), ...
                    'LineWidth', 2)];
            elseif i_fig==2
                a2 = [a2; plot3(y(interval,i_fig),y(interval,1+i_fig),y(interval,2+i_fig), ...
                    'LineWidth', 2)];
            end
            hold on
            
            
            scatter3(y(interval,i_fig),y(interval,1+i_fig),y(interval,2+i_fig), ...
                15,c1,'filled')
            

            xlabel(['PC',num2str(i_fig)])
            ylabel(['PC',num2str(1+i_fig)])
            zlabel(['PC',num2str(2+i_fig)])
            title(stimtypes{i_stim})
            
            
        end
    end
    
    figure(3646*1+i_stim); legend(a1, a1a, 'Location','NorthEast')
    figure(3646*2+i_stim); legend(a2, a1a, 'Location','NorthEast')

end

%% save figs

pathout = [FileIn_path, '\', Series_ID, '\figures\', mfilename, ...
    '\from', num2str(intv_start), 'to', num2str(intv_end), 's'];
if ~exist(pathout,'dir'); mkdir(pathout); end

nameflag = 'alltrials';
FileOut = [pathout,'\', nameflag,'_PC123.fig']
savefig(figure(2865),FileOut)
FileOut = [pathout,'\', nameflag,'_PC234.fig']
savefig(figure(2866),FileOut)


for i_stim = 1:nstimtypes
    
    nameflag = [stimtypes{i_stim}];
    
    FileOut = [pathout,'\', nameflag,'_PC123.fig']
    savefig(figure(3646*1+i_stim),FileOut)
    
    FileOut = [pathout,'\', nameflag,'_PC234.fig']
    savefig(figure(3646*2+i_stim),FileOut)
end




