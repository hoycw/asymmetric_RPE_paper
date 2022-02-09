function SBJ11e_HFA_conn_plot_grp_mLME_confusion(proc_id, an_id, model_id, conn_id, swap_Xcorr, stat_id)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
[reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
[~, names, colors, ~, ~] = fn_puns_category_label_styles('puns');

%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/conn/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir proc_id '_' model_id '_' an_id '_' conn_id '_chancoef.mat'];
load(stats_fname)

if swap_Xcorr == 1
    for r = 1:numel(conn_stats_chan.coefs)
        conn_stats_chan.pair_label{r} = conn_stats_chan.pair_label{r}(1,[2,1]);
        for b = 1:size(conn_stats_chan.coefs{r},4)
            for ch = 1:size(conn_stats_chan.coefs{r},1)
                conn_stats_chan.coefs{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.coefs{r}(ch,:,:,b)));
                conn_stats_chan.pvals{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.pvals{r}(ch,:,:,b)));
                conn_stats_chan.qvals{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.qvals{r}(ch,:,:,b)));
                conn_stats_chan.lower{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.lower{r}(ch,:,:,b)));
                conn_stats_chan.upper{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.upper{r}(ch,:,:,b)));
                %conn_stats_chan.chan_label{r} = conn_stats_chan.chan_label{r}([2,1]);
            end
        end
    end
end

%%
hfa_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(hfa_fname)
%% select and classify channel pairs
confarray = cell(numel(beta_chan.chancat_label));
maskarray = cell(numel(beta_chan.chancat_label));
for b = 1:size(beta_chan.coefs{1},4)
    for cat1 = 1:numel(beta_chan.chancat_label)
        for cat2 = 1:numel(beta_chan.chancat_label)
            for ch1 = 1:length(beta_chan.chancat_ix{1}{cat1})
                for ch2 = 1:length(beta_chan.chancat_ix{2}{cat2})
                    ch_lab1 = beta_chan.chan_label{1}{beta_chan.chancat_ix{1}{cat1}(ch1)};
                    ch_lab2 = beta_chan.chan_label{2}{beta_chan.chancat_ix{2}{cat2}(ch2)};
                    if strcmp(ch_lab1(1:4), ch_lab2(1:4))
                        chan_pair = [ch_lab1 '_to_' ch_lab2(6:end)];
                        display(chan_pair)
                        conn_ix = find(strcmp(chan_pair, conn_stats_chan.chan_label{1}));
                        display(conn_stats_chan.chan_label{1}{conn_ix})
                        if ~isempty(conn_ix)
                            confarray{cat1,cat2}(end+1,:,:) = conn_stats_chan.coefs{1}(conn_ix,:,:,1);
                            maskarray{cat1,cat2}(end+1,:,:) = conn_stats_chan.qvals{1}(conn_ix,:,:,1);
                        end
                    end
                end
            end
        end
    end
end
%% plot confusion matrix
for rg = 3:4
    plotix = 0;
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    for cr = 1:size(confarray,1)
        for cc = 1:size(confarray,2)
            plotix = plotix +1;
            subplot(size(confarray,1),size(confarray,2), plotix)
            yline(0); hold on; xline(0);hold on;
            ylim([-0.07,0.07])
            if ~isempty(confarray{cr,cc})
                qmask = double(squeeze(maskarray{cr,cc}(:,rg,:) < 0.05));
                qmask(qmask == 0) = NaN;
                plot(conn_stats_chan.time, squeeze(confarray{cr,cc}(:,rg,:)),'k');hold on;
                plot(conn_stats_chan.time, squeeze(confarray{cr,cc}(:,rg,:)).*qmask,'k',...
                'LineWidth',1.5);

                if cr == 1
                    title([beta_chan.chancat_label{cc} ' ' beta_chan.label{2}])
                end
                if cc == 1
                    ylabel([beta_chan.chancat_label{cr} ' ' beta_chan.label{1}])
                end
            end
        end
    end
    sgtitle([reg_lab{rg-1} ' coefficients'])
    plot_fname = sprintf('%s%s_%s_%s_hfa_confusion_%s_to_%s_%s.pdf', fig_dir,...
                proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
                conn_stats_chan.pair_label{r}{2},reg_lab{rg-1});
        print(plot_fname,cf,'-dpdf','-fillpage')
end
%% Confusion plot
% ptime = 0.05;
% t = find(round(conn_stats_chan.time,3) == ptime);
% confmeans = NaN(size(confarray));
% for cr = 1:size(confarray,1)
%     for cc = 1:size(confarray,2)
%         if ~isempty(confarray{cr,cc})
%             confmeans(cr,cc) = squeeze(mean(abs(confarray{cr,cc}(:,4,t)),1));
%         end
%     end
% end
% figure;
% imagesc(confmeans)
% colorbar()

% ntimes = numel(conn_stats_chan.time);
% for r = 1:numel(conn_stats_chan.coefs)
%     for b = 1:size(conn_stats_chan.coefs{r},4)
%         cf = figure('units','normalized','outerposition',[0 0 1 1],...
%             'PaperOrientation','Landscape');
%         pchan = conn_stats_chan.chancat_ix{r}{1,b};
%         nchan = conn_stats_chan.chancat_ix{r}{2,b};
%         slnchan = conn_stats_chan.chancat_ix{r}{3,b};
%         rwdchan = conn_stats_chan.chancat_ix{r}{4,b};
%         for t = 1:ntimes
%             cplots = [];
%             subplot(ceil(ntimes / 5),5,t)
%             yline(0); hold on; xline(0);
%             for ctg = 1:numel(conn_stats_chan.chancat_label)
%                 cp = scatter(conn_stats_chan.coefs{r}(conn_stats_chan.chancat_ix{r}{ctg,b},3,t,b),...
%                     conn_stats_chan.coefs{r}(conn_stats_chan.chancat_ix{r}{ctg,b},4,t,b),...
%                     [],colors{ctg},'filled','p','LineWidth',1,...
%                     'MarkerEdgeColor',colors{ctg});
%                 cp.MarkerFaceAlpha = 0.3;
%                 cp.MarkerEdgeAlpha = 0.7;
%                 cplots(ctg) =  cp;
%                 hold on;
%             end
%             
%             title(sprintf('%.00f ms lag',conn_stats_chan.time(t)*1000));
%             
%             xlim([-0.06,0.06])
%             ylim([-0.04,0.04])
%             
%             if t == length(conn_stats_chan.time)
%                 xlabel('pRPE coefficient (a.u.)', 'FontSize',7)
%                 ylabel('nRPE coefficient (a.u.)', 'FontSize',7)
%                 hold off;
%                 
%                 legend(cplots,names,'FontSize',9,'Location','northeast',...
%                     'NumColumns',4, 'Position',[0.45,0.02,0.1,0.03])
%                 %legend('boxoff')
%             end
%         end
%         sgtitle([conn_stats_chan.pair_label{r}{1} ' to ' conn_stats_chan.pair_label{r}{2}])
%         plot_fname = sprintf('%s%s_%s_%s_hfa_chanscatter_%s_to_%s_%d.pdf', fig_dir,...
%             proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
%             conn_stats_chan.pair_label{r}{2},b);
%         print(plot_fname,cf,'-dpdf','-fillpage')
%     end
% end
% %close all
% %% plot static scatterplot
% ptime = 0.05;
% t = find(round(conn_stats_chan.time,3) == ptime);
% ntimes = numel(conn_stats_chan.time);
% 
% for r = 1:numel(conn_stats_chan.coefs)
%     for b = 1:size(conn_stats_chan.coefs{r},4)
%         cf = figure('units','normalized','outerposition',[0 0 1 1],...
%             'PaperOrientation','Landscape');
%         pchan = conn_stats_chan.chancat_ix{r}{1,b};
%         nchan = conn_stats_chan.chancat_ix{r}{2,b};
%         slnchan = conn_stats_chan.chancat_ix{r}{3,b};
%         rwdchan = conn_stats_chan.chancat_ix{r}{4,b};
%         cplots = [];
%         yline(0); hold on; xline(0);
%         for ctg = 1:numel(conn_stats_chan.chancat_label)
%             cp = scatter(conn_stats_chan.coefs{r}(conn_stats_chan.chancat_ix{r}{ctg,b},3,t,b),...
%                 conn_stats_chan.coefs{r}(conn_stats_chan.chancat_ix{r}{ctg,1},4,t,b),...
%                 150,colors{ctg},'filled','p','LineWidth',1,...
%                 'MarkerEdgeColor',colors{ctg});
%             cp.MarkerFaceAlpha = 0.3;
%             cp.MarkerEdgeAlpha = 1;
%             
%             cplots(ctg) =  cp;
%             hold on;
%         end
%         
%         title(sprintf('%s to %s (%.00f ms lag)',conn_stats_chan.pair_label{r}{1},...
%             conn_stats_chan.pair_label{r}{2},conn_stats_chan.time(t)*1000));
%         xlim([-0.06,0.06])
%         ylim([-0.04,0.04])
%         
%         
%         xlabel('Positive RPE coefficient (a.u.)', 'FontSize',11)
%         ylabel('Negative RPE coefficient (a.u.)', 'FontSize',11)
%         hold off;
%         
%         %legend('boxoff')
%         legend(cplots,names,'FontSize',9,...%'Location','northeast',...
%             'NumColumns',1, 'Position',[0.15,0.7,0.1,0.1])
%         plot_fname = sprintf('%s%s_%s_%s_hfa_chanscatter_%s_to_%s_%d_static.pdf', fig_dir,...
%             proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
%             conn_stats_chan.pair_label{r}{2},b);
%         print(plot_fname,cf,'-dpdf','-fillpage')
%        
%     end
% end
% close all
% %% Plot channel predicted values 
% % coef_lims = [-1.7,1.7];
% % for r = 1:numel(conn_stats_chan.coefs)
% %     for b = 1:size(conn_stats_chan.coefs{r},5)
% %         for cl = 1:length(coef_lims)
% %             cf = figure('units','normalized','outerposition',[0 0 1 1],...
% %                 'PaperOrientation','Landscape');
% %             for rg = 1:length(reg_lab)
% %                 tcourses1 = squeeze(conn_stats_chan.coefs{r}(:,1,:,b));% + conn_stats_chan.coefs{r}(:,rg+1,:,b) * -1.7;
% %                 tcourses2 = squeeze(conn_stats_chan.coefs{r}(:,rg+1,:,b)) * coef_lims(cl);
% %                 [cchans,~] = find(tcourses2 < 0 & (tcourses1 + tcourses2) < 0);
% %                 cchans = unique(cchans);
% %                 subplot(1,3,rg)
% %                 yline(0); hold on; xline(0); hold on;
% %                 hga1 = plot(conn_stats_chan.time, tcourses1(cchans,:) + tcourses2(cchans,:), 'color','k');hold on;
% %                 %hga2 = plot(conn_stats_chan.time, tcourses1(cchans,:), 'color','k');
% %                 for hg = 1:length(hga1)
% %                     hga1(hg).Color(4) = 0.1;  hold on;
% %                     %hga2(hg).Color(4) = 0.1;  hold on;
% %                 end
% %                 ylim([-0.1,0.1])
% %             end
% %         end
% %     end
% % end
% %% Plot channel predicted values (to be continued)
% ptime = 0.05;
% t = find(round(conn_stats_chan.time,3) == ptime);
% for r = 1:numel(conn_stats_chan.coefs)
%     for b = 1:size(conn_stats_chan.coefs{r},4)
%         tcourses1 = squeeze(conn_stats_chan.coefs{r}(:,1,t,b));% + conn_stats_chan.coefs{r}(:,rg+1,:,b) * -1.7;
%         cf = figure('units','normalized','outerposition',[0 0 1 1],...
%             'PaperOrientation','Landscape');
%         for rg = 2:length(reg_lab)
%             tcourses2 = squeeze(conn_stats_chan.coefs{r}(:,rg+1,t,1)) * 1.7;
%             %nidx = find((tcourses1 + tcourses2) < 0 & tcourses2 < 0 & conn_stats_chan.qvals{1}(:,rg+1,t,1) < .05);
%             nidx = find(conn_stats_chan.qvals{r}(:,rg+1,t,1) < .05); % significant channels
%             subplot(1,2,rg-1)
%             yline(0); hold on;
%             hga1 = plot([0,1],[tcourses1, tcourses1 + tcourses2],'.k');hold on;
%             hga2 = plot([0,1],[tcourses1, tcourses1 + tcourses2],'-k');hold on;
%             
%             for hg = 1:length(hga1)
%                 hga1(hg).Color(4) = 0.05;
%                 hga2(hg).Color(4) = 0.05;
%             end
%             if ~isempty(nidx)
%                 hga3 = plot([0,1],[tcourses1(nidx), tcourses1(nidx) + tcourses2(nidx)],'.r');hold on;
%                 hga4 = plot([0,1],[tcourses1(nidx), tcourses1(nidx) + tcourses2(nidx)],'-r');hold on;
%                 
%                 for hg = 1:length(hga3)
%                     hga3(hg).Color(4) = 0.5;
%                     hga4(hg).Color(4) = 0.5;
%                 end
%             end
%             xlim([-0.5,1.5])
%             ylim([-0.15,0.25])
%             ax = gca;
%             ax.XTick = [0,1];
%             ax.XTickLabels = {'Intercept',['Intercept + ' reg_lab{rg} ' * 1.7' ]};
%             ylabel('cross-correlation coefficient (r)')
%             if rg == 3
%             legend([hga1(1),hga3(1)],{'all channel pairs','significant channel pairs'},...
%                     'FontSize',9,'Location','northeast',...
%                     'NumColumns',1)%, 'Position',[0.45,0.02,0.1,0.03])
%             end
%         end
%         sgtitle(sprintf('Predicted connectivity from %s to %s at %.00f ms lag',...
%                         conn_stats_chan.pair_label{r}{1},...
%                         conn_stats_chan.pair_label{r}{2}, ptime*1000))
%         plot_fname = sprintf('%s%s_%s_%s_hfa_chancoef_%s_to_%s_%d_prediction.pdf', fig_dir,...
%                         proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
%                          conn_stats_chan.pair_label{r}{2},b);
%         print(plot_fname,cf,'-dpdf','-fillpage')
%     end
% end
% %% plot channel pairs time courses
% for r = 1:numel(conn_stats_chan.coefs)
%     % get mask
%     for b = 1:size(conn_stats_chan.coefs{r},5)
%         cf = figure('units','normalized','outerposition',[0 0 1 1],...
%         'PaperOrientation','Landscape');
%         cplots = NaN(1,length(conn_stats_chan.chancat_label));
%         for rg = 1:length(reg_lab)
%             [ridx,~] = find(squeeze(conn_stats_chan.qvals{r}(:, rg + 1,:,b) < .05));
%             ylims0 = [min(conn_stats_chan.coefs{r}(:,2:end,:,b),[],'all'),...
%                 max(conn_stats_chan.coefs{r}(:,2:end,:,b),[],'all')];
%             ylims0(1) = min([ylims0(1) - 0.2*abs(ylims0(1)), 0]);
%             ylims0(2) = max([ylims0(2) + 0.2*abs(ylims0(2)), 0]);
%             ylims = [min(ylims0(1)),max(ylims0(2))];
%             subplot(1,3,rg)
%             yline(0); hold on; xline(0); hold on;
%             for ch = 1:size(conn_stats_chan.coefs{r},1)
%                 qmask = squeeze(double(conn_stats_chan.pvals{r}(ch, rg + 1,:) < .05));
%                 qmask(qmask == 0) = Inf;
%                 ch_cat = [];
%                 for ctg = 1:length(conn_stats_chan.chancat_label)
%                     ch_cat(ctg) = ismember(ch, conn_stats_chan.chancat_ix{r}{ctg});
%                 end
%                 if sum(ch_cat) > 0
%                     if group_colors == 1
%                         ccolor = colors{find(ch_cat)};
%                     else
%                         ccolor = 'k';
%                     end
%                     hga = plot(conn_stats_chan.time, squeeze(conn_stats_chan.coefs{r}(ch,rg+1,:,b)),...
%                         'color', ccolor); hold on;
%                     hga2 = plot(conn_stats_chan.time, squeeze(conn_stats_chan.coefs{r}(ch,rg+1,:,b)).*qmask,...
%                         'color',ccolor,'LineWidth',1.5); hold on;
%                     cplots(find(ch_cat)) = hga2;
%                 else
%                     hga = plot(conn_stats_chan.time, squeeze(conn_stats_chan.coefs{r}(ch, rg + 1,:,b)),'k-');
%                     hga.Color(4) = 0.1;  hold on;
%                 end
%                 title(reg_names{rg})
%                 ylim(ylims)
%                 xlabel('time lag (s)')
%                 ylabel('coefficient (a.u.)')
%             end
%         end
%         if group_colors == 1
%             legend(cplots,names,'FontSize',9,... %'Location','northeast',...
%                 'NumColumns',1, 'Position',[0.8,0.75,0.05,0.1])
%         end
%         sgtitle([conn_stats_chan.pair_label{r}{1} ' to ' conn_stats_chan.pair_label{r}{2}])
%         plot_fname = sprintf('%s%s_%s_%s_hfa_chancoef_%s_to_%s_%d.pdf', fig_dir,...
%             proc_id, model_id, an_id, conn_stats_chan.pair_label{r}{1},...
%             conn_stats_chan.pair_label{r}{2},b);
%         print(plot_fname,cf,'-dpdf','-fillpage')
%     end
% end
% close all