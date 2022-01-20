function SBJ08h_HFA_plot_grp_mLME_chancoef(proc_id, an_id, model_id, stat_id)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
[reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname)

%% plot channel timecourses
for r = 1:numel(beta_chan.coefs)
    % get mask
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    
    for rg = 1:length(reg_lab)
        [ridx,~] = find(squeeze(beta_chan.qvals{r}(:, rg + 1,:) < .05));
        ylims0 = [min(beta_chan.coefs{r}(:,r+1,:),[],'all'),...
            max(beta_chan.coefs{r}(:,r+1,:),[],'all')];
        ylims0(1) = min([ylims0(1) - 0.5*abs(ylims0(1)), 0]);
        ylims0(2) = max([ylims0(2) + 0.5*abs(ylims0(2)), 0]);
        ylims = [min(ylims0(1)),max(ylims0(2))];
        subplot(1,3,rg)
        yline(0); hold on; xline(0); hold on;
        for ch = 1:size(beta_chan.coefs{r},1)
            %qmask = squeeze(double(beta_chan.pvals{r}(ch, rg + 1,:) < .05));
            %qmask(qmask == 0) = Inf;
            hga = plot(beta_chan.time, squeeze(beta_chan.coefs{r}(ch,rg+1,:)),'k-'); hold on;
            hga.Color(4) = 0.1;  hold on;
            if ismember(ch,ridx)
                plot(beta_chan.time, squeeze(beta_chan.coefs{r}(ch, rg + 1,:)),'k-')
            end
        end
        title(reg_names{rg})
        ylim(ylims)
        xlabel('time (s)')
        ylabel('coefficient (a.u.)')
    end
    sgtitle(beta_chan.label{r})
    plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa_chancoef_' beta_chan.label{r} '.pdf'];
    print(plot_fname,cf,'-dpdf','-fillpage')
end
close all


%% plot dynamic scatterplot
alphav = 0.3;
ntimes = numel(beta_chan.time);
for r = 1:numel(beta_chan.coefs)
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    [pidx,~] = find(squeeze(beta_chan.qvals{r}(:,3,:) < .05));
    [nidx,~] = find(squeeze(beta_chan.qvals{r}(:,4,:) < .05));
    pidx = unique(pidx); nidx = unique(nidx);
    pchan = pidx(~ismember(pidx,nidx));
    nchan = nidx(~ismember(nidx,pidx));
    for t = 1:ntimes
        signcoef = double(beta_chan.coefs{r}(:,3:4,t) >= 0);
        slnidx = find(signcoef(:,1) ~= signcoef(:,2));
        rwdidx = find(signcoef(:,1) == signcoef(:,2));
        slnchan = pidx(ismember(pidx,nidx));
        slnchan = slnchan(ismember(slnchan,slnidx));
        rwdchan = pidx(ismember(pidx,nidx));
        rwdchan = rwdchan(ismember(rwdchan,rwdidx));
        
        subplot(ceil(ntimes / 5),5,t)
        s1= scatter(beta_chan.coefs{r}(pchan,3,t),beta_chan.coefs{r}(pchan,4,t)*-1,[],...
            'g','filled', 'LineWidth',1,'MarkerEdgeColor','k');
        s1.MarkerFaceAlpha = alphav; s1.MarkerEdgeAlpha = alphav; hold on;
        s2= scatter(beta_chan.coefs{r}(nchan,3,t),beta_chan.coefs{r}(nchan,4,t)*-1,[],...
            'y','filled', 'LineWidth',1,'MarkerEdgeColor','k');hold on;
        s2.MarkerFaceAlpha = alphav; s2.MarkerEdgeAlpha = alphav; hold on;
        s3= scatter(beta_chan.coefs{r}(slnchan,3,t),beta_chan.coefs{r}(slnchan,4,t)*-1,[],...
            'r','filled', 'LineWidth',1,'MarkerEdgeColor','k');hold on;
        s3.MarkerFaceAlpha = alphav; s3.MarkerEdgeAlpha = alphav; hold on;
        s4= scatter(beta_chan.coefs{r}(rwdchan,3,t),beta_chan.coefs{r}(rwdchan,4,t)*-1,[],...
            'b','filled', 'LineWidth',1,'MarkerEdgeColor','k');
        s4.MarkerFaceAlpha = alphav; s4.MarkerEdgeAlpha = alphav; hold on;
        
        yline(0); hold on; xline(0);
        title(sprintf('%.00f ms',beta_chan.time(t)*1000));
        xlim([-0.7,0.7])
        ylim([-1,1])
        
        if t == length(beta_chan.time)
            xlabel('pRPE coefficient (a.u.)', 'FontSize',7)
            ylabel('(-) nRPE coefficient (a.u.)', 'FontSize',7)
            hold off;
            
            legend([s1,s2,s3,s4],{'pRPE','nRPE','uRPE','sRPE'},'FontSize',9,...%'Location','northeast',...
                'NumColumns',1, 'Position',[0.625,0.1,0.1,0.1])
            %legend('boxoff')
        end
    end
    sgtitle(beta_chan.label{r})
    plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa_chanscatter_' beta_chan.label{r} '.pdf'];
    print(plot_fname,cf,'-dpdf','-fillpage')
end
close all
