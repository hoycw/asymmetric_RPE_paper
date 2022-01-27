function SBJ11d_HFA_conn_plot_grp_mLME_chancoef(proc_id, an_id, model_id, conn_id)

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
fig_dir = [root_dir 'PRJ_Error/results/conn/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir proc_id '_' model_id '_' an_id '_' conn_id '_chancoef.mat'];
load(stats_fname)

%% plot channel pairs time courses
for r = 1:numel(conn_stats_chan.coefs)
    % get mask
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    for b = 1:size(conn_stats_chan.coefs{r},4)
        for rg = 1:length(reg_lab)
            [ridx,~] = find(squeeze(conn_stats_chan.qvals{r}(:, rg + 1,:,b) < .05));
            ylims0 = [min(conn_stats_chan.coefs{r}(:,r+1,:,b),[],'all'),...
                max(conn_stats_chan.coefs{r}(:,r+1,:,b),[],'all')];
            ylims0(1) = min([ylims0(1) - 0.8*abs(ylims0(1)), 0]);
            ylims0(2) = max([ylims0(2) + 0.8*abs(ylims0(2)), 0]);
            ylims = [min(ylims0(1)),max(ylims0(2))];
            subplot(1,3,rg)
            yline(0); hold on; xline(0); hold on;
            for ch = 1:size(conn_stats_chan.coefs{r},1)
                %qmask = squeeze(double(beta_chan.pvals{r}(ch, rg + 1,:) < .05));
                %qmask(qmask == 0) = Inf;
                hga = plot(conn_stats_chan.time, squeeze(conn_stats_chan.coefs{r}(ch,rg+1,:,b)),'k-'); hold on;
                hga.Color(4) = 0.05;  hold on;
                if ismember(ch,ridx)
                    plot(conn_stats_chan.time, squeeze(conn_stats_chan.coefs{r}(ch, rg + 1,:,b)),'k-')
                end
            end
            title(reg_names{rg})
            ylim(ylims)
            xlabel('time lag (s)')
            ylabel('coefficient (a.u.)')
        end
        sgtitle([conn_stats_chan.pair_label{r}{1} ' to ' conn_stats_chan.pair_label{r}{2}])
        plot_fname = sprintf('%s%s_%s_%s_hfa_chancoef_%s_to_%s_%d.pdf', fig_dir,...
            proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
            conn_stats_chan.pair_label{r}{2},b);
        print(plot_fname,cf,'-dpdf','-fillpage')
    end
end
%close all


%% plot dynamic scatterplot
alphav = 0.3;
ntimes = numel(conn_stats_chan.time);
for r = 1:numel(conn_stats_chan.coefs)
    for b = 1:size(conn_stats_chan.coefs{r},4)
        cf = figure('units','normalized','outerposition',[0 0 1 1],...
            'PaperOrientation','Landscape');
        [pidx,~] = find(squeeze(conn_stats_chan.qvals{r}(:,3,:,b) < .05));
        [nidx,~] = find(squeeze(conn_stats_chan.qvals{r}(:,4,:,b) < .05));
        pidx = unique(pidx); nidx = unique(nidx);
        pchan = pidx(~ismember(pidx,nidx));
        nchan = nidx(~ismember(nidx,pidx));
        for t = 1:ntimes
            signcoef = double(conn_stats_chan.coefs{r}(:,3:4,t,b) >= 0);
            slnidx = find(signcoef(:,1) ~= signcoef(:,2));
            rwdidx = find(signcoef(:,1) == signcoef(:,2));
            slnchan = pidx(ismember(pidx,nidx));
            slnchan = slnchan(ismember(slnchan,slnidx));
            rwdchan = pidx(ismember(pidx,nidx));
            rwdchan = rwdchan(ismember(rwdchan,rwdidx));
            
            subplot(ceil(ntimes / 5),5,t)
            s1= scatter(conn_stats_chan.coefs{r}(pchan,3,t,b),conn_stats_chan.coefs{r}(pchan,4,t,b)*-1,[],...
                'g','filled', 'LineWidth',1,'MarkerEdgeColor','k');
            s1.MarkerFaceAlpha = alphav; s1.MarkerEdgeAlpha = alphav; hold on;
            s2= scatter(conn_stats_chan.coefs{r}(nchan,3,t,b),conn_stats_chan.coefs{r}(nchan,4,t,b)*-1,[],...
                'y','filled', 'LineWidth',1,'MarkerEdgeColor','k');hold on;
            s2.MarkerFaceAlpha = alphav; s2.MarkerEdgeAlpha = alphav; hold on;
            s3= scatter(conn_stats_chan.coefs{r}(slnchan,3,t,b),conn_stats_chan.coefs{r}(slnchan,4,t,b)*-1,[],...
                'r','filled', 'LineWidth',1,'MarkerEdgeColor','k');hold on;
            s3.MarkerFaceAlpha = alphav; s3.MarkerEdgeAlpha = alphav; hold on;
            s4= scatter(conn_stats_chan.coefs{r}(rwdchan,3,t,b),conn_stats_chan.coefs{r}(rwdchan,4,t,b)*-1,[],...
                'b','filled', 'LineWidth',1,'MarkerEdgeColor','k');
            s4.MarkerFaceAlpha = alphav; s4.MarkerEdgeAlpha = alphav; hold on;
            
            yline(0); hold on; xline(0);
            title(sprintf('%.00f ms lag',conn_stats_chan.time(t)*1000));
            xlim([-0.06,0.06])
            ylim([-0.04,0.04])
            
            if t == length(conn_stats_chan.time)
                xlabel('pRPE coefficient (a.u.)', 'FontSize',7)
                ylabel('(-) nRPE coefficient (a.u.)', 'FontSize',7)
                hold off;
                
                legend([s1,s2,s3,s4],{'pRPE','nRPE','uRPE','sRPE'},'FontSize',9,...%'Location','northeast',...
                    'NumColumns',4, 'Position',[0.45,0.02,0.1,0.03])
                %legend('boxoff')
            end
        end
        sgtitle([conn_stats_chan.pair_label{r}{1} ' to ' conn_stats_chan.pair_label{r}{2}])
        plot_fname = sprintf('%s%s_%s_%s_hfa_chanscatter_%s_to_%s_%d.pdf', fig_dir,...
            proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
            conn_stats_chan.pair_label{r}{2},b);
        print(plot_fname,cf,'-dpdf','-fillpage')

    end
end
close all
