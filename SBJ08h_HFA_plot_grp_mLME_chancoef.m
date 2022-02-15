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

if isfield(beta_chan,'chancat_label')
    [~, names, colors, ~, ~] = fn_puns_category_label_styles('puns');
else
    names = {'significant channels'};
    colors = {'k'};
    for r = 1:numel(beta_chan.coefs)
        beta_chan.chancat_label = {'sigchan'};
        [sig_ix, ~] = find(beta_chan.qvals{r}(:,3:end,:) < 0.05);
        beta_chan.chancat_ix{r} = {unique(sig_ix)};
    end
end
%% plot channel time courses
for r = 1:numel(beta_chan.coefs)
    % get mask
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
                'PaperOrientation','Landscape');
                %'PaperOrientation','Landscape')
    cplots = NaN(1,length(beta_chan.chancat_label));
    for rg = 1:length(reg_lab)
        %[ridx,~] = find(squeeze(beta_chan.qvals{r}(:, rg + 1,:) < .05));
        ylims0 = [min(beta_chan.coefs{r}(:,2:end,:),[],'all'),...
            max(beta_chan.coefs{r}(:,2:end,:),[],'all')];
        ylims0(1) = min([ylims0(1) - 0.5*abs(ylims0(1)), 0]);
        ylims0(2) = max([ylims0(2) + 0.5*abs(ylims0(2)), 0]);
        ylims = [min(ylims0(1)),max(ylims0(2))];
        subplot(1,3,rg)
        yline(0); hold on; xline(0); hold on;
        for ch = 1:size(beta_chan.coefs{r},1)
            qmask = squeeze(double(beta_chan.qvals{r}(ch, rg + 1,:) < .05));
            qmask(qmask == 0) = Inf;
            ch_cat = [];
            for ctg = 1:length(beta_chan.chancat_label)
                ch_cat(ctg) = ismember(ch, beta_chan.chancat_ix{r}{ctg});
            end
            if sum(ch_cat) > 0
                hga = plot(beta_chan.time, squeeze(beta_chan.coefs{r}(ch,rg+1,:)),...
                    'color', colors{find(ch_cat)}); hold on;
                hga2 = plot(beta_chan.time, squeeze(beta_chan.coefs{r}(ch,rg+1,:)).*qmask,...
                       'color',colors{find(ch_cat)},'LineWidth',1.5); hold on;
                hga3 = plot(beta_chan.time, squeeze(beta_chan.coefs{r}(ch,rg+1,:)).*qmask,...
                       '.','color',colors{find(ch_cat)},'LineWidth',1.25); hold on;
                cplots(find(ch_cat)) = hga2;
            else
                hga = plot(beta_chan.time, squeeze(beta_chan.coefs{r}(ch, rg + 1,:)),'k-');
                hga.Color(4) = 0.1;  hold on;
            end
        end
        title(reg_names{rg})
        ylim(ylims)
        xlabel('time (s)')
        ylabel('coefficient (a.u.)')
    end
    legend(cplots,names,'FontSize',9,...%'Location','northeast',...
        'NumColumns',1, 'Position',[0.8,0.75,0.05,0.1])
    sgtitle(beta_chan.label{r})
    plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa_chancoef_' beta_chan.label{r} '.pdf'];
    print(plot_fname,cf,'-dpdf','-fillpage')
end
close all

if strcmp(model_id,'EpnRPE_DifFB')
    %% plot dynamic scatterplot
    ntimes = numel(beta_chan.time);
    for r = 1:numel(beta_chan.coefs)
        cf = figure('units','normalized','outerposition',[0 0 1 1],...
            'PaperOrientation','Landscape');
        pchan = beta_chan.chancat_ix{r}{1,1};
        nchan = beta_chan.chancat_ix{r}{2,1};
        slnchan = beta_chan.chancat_ix{r}{3,1};
        rwdchan = beta_chan.chancat_ix{r}{4,1};
        for t = 1:ntimes
            cplots = [];
            subplot(ceil(ntimes / 5),5,t)
            yline(0); hold on; xline(0);
            for ctg = 1:numel(beta_chan.chancat_label)
                cp = scatter(beta_chan.coefs{r}(beta_chan.chancat_ix{r}{ctg,1},3,t),...
                    beta_chan.coefs{r}(beta_chan.chancat_ix{r}{ctg,1},4,t),...
                    [],colors{ctg},'filled','p','LineWidth',1,...
                    'MarkerEdgeColor',colors{ctg});
                cp.MarkerFaceAlpha = 0.3;
                cp.MarkerEdgeAlpha = 0.7;
                cplots(ctg) =  cp;
                hold on;
            end
            
            title(sprintf('%.00f ms',beta_chan.time(t)*1000));
            xlim([-0.7,0.7])
            ylim([-0.8,0.8])
            
            if t == length(beta_chan.time)
                xlabel('pRPE coefficient (a.u.)', 'FontSize',7)
                ylabel('(-) nRPE coefficient (a.u.)', 'FontSize',7)
                hold off;
                
                legend(cplots,names,'FontSize',9,...%'Location','northeast',...
                    'NumColumns',1, 'Position',[0.625,0.1,0.1,0.1])
                %legend('boxoff')
            end
        end
        sgtitle(beta_chan.label{r})
        plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa_chanscatter_' beta_chan.label{r} '.pdf'];
        print(plot_fname,cf,'-dpdf','-fillpage')
    end
    close all
    %% plot static scatterplot
    ptime = 0.350;
    t = find(round(beta_chan.time,3) == ptime);
    ntimes = numel(beta_chan.time);
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    for r = 1:numel(beta_chan.coefs)
        pchan = beta_chan.chancat_ix{r}{1,1};
        nchan = beta_chan.chancat_ix{r}{2,1};
        slnchan = beta_chan.chancat_ix{r}{3,1};
        rwdchan = beta_chan.chancat_ix{r}{4,1};
        subplot(1,2,r)
        cplots = [];
        yline(0); hold on; xline(0);
        for ctg = 1:numel(beta_chan.chancat_label)
            cp = scatter(beta_chan.coefs{r}(beta_chan.chancat_ix{r}{ctg,1},3,t),...
                beta_chan.coefs{r}(beta_chan.chancat_ix{r}{ctg,1},4,t),...
                150,colors{ctg},'filled','p','LineWidth',1,...
                'MarkerEdgeColor',colors{ctg});
            cp.MarkerFaceAlpha = 0.3;
            cp.MarkerEdgeAlpha = 1;
            
            cplots(ctg) =  cp;
            hold on;
        end
        
        title(sprintf('%s (%.00f ms)',beta_chan.label{r},beta_chan.time(t)*1000));
        xlim([-0.7,0.7])
        ylim([-0.8,0.8])
        
        
        xlabel('Positive RPE coefficient (a.u.)', 'FontSize',11)
        ylabel('Negative RPE coefficient (a.u.)', 'FontSize',11)
        hold off;
        
        %legend('boxoff')
        
    end
    legend(cplots,names,'FontSize',9,...%'Location','northeast',...
        'NumColumns',1, 'Position',[0.15,0.7,0.1,0.1])
    plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa_chanscatter_' beta_chan.label{r} '_static.pdf'];
    print(plot_fname,cf,'-dpdf','-fillpage')
    close all
    
    %% Plot trajectories
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    ntimes = numel(beta_chan.time);
    for r = 1:numel(beta_chan.coefs)
        
        pchan = beta_chan.chancat_ix{r}{1,1};
        nchan = beta_chan.chancat_ix{r}{2,1};
        slnchan = beta_chan.chancat_ix{r}{3,1};
        rwdchan = beta_chan.chancat_ix{r}{4,1};
        
        cplots = [];
        subplot(1,2,r)
        yline(0); hold on; xline(0);
        for ctg = 1:numel(beta_chan.chancat_label)
            time3d = linspace(beta_chan.time(1),beta_chan.time(end));
            x = interp1(beta_chan.time, squeeze(beta_chan.coefs{r}(beta_chan.chancat_ix{r}{ctg,1},3,:))',...
                time3d,'spline');
            y = interp1(beta_chan.time, squeeze(beta_chan.coefs{r}(beta_chan.chancat_ix{r}{ctg,1},4,:))',...
                time3d,'spline');
            if size(x,2) == 100
                x = x';
                y = y';
            end
            cp = plot(x,y,'color',colors{ctg},'LineWidth',1); hold on;
            %         cp2 = scatter(x(1,:),y(1,:),35,colors{ctg},'filled','o',...
            %         'MarkerEdgeColor','k');
            %         cp2.MarkerFaceAlpha = 0.8;
            %         cp2.MarkerEdgeAlpha = 0.8;
            %        hold on;
            for nc = 1:size(x,2)
                arrtimes = [1,40,50,60,70,99];
                for ar = 1:length(arrtimes)
                    cart = arrtimes(ar);
                    ah = annotation('arrow','color',colors{ctg},...
                        'headStyle','cback1','HeadLength',3,'HeadWidth',5);
                    set(ah,'parent',gca);
                    set(ah,'position',[x(cart,nc), y(cart,nc), x(cart+1,nc) - x(cart,nc), y(cart+1,nc) - y(cart,nc)]);
                end
            end
            %         quiver(x(50,:), y(50,:), x(51,:) - x(50,:) , y(51,:) - y(50,:),...
            %                'off','color',colors{ctg},'LineStyle','none','ShowArrowHead','on')
            %cp = plot3(x,time3d,y,'color',colors{ctg},'LineWidth',1.25); grid on;
            %         for cpp = 1:length(cp)
            %             cp(cpp).Color(4) = 0.6;
            %         end
            cplots(ctg) =  cp(1);
            hold on;
        end
        
        xlim([-0.35,0.6])
        ylim([-0.35,0.7])
        %ylim([time3d(1),time3d(end)])
        xlabel('pRPE coefficient (a.u.)', 'FontSize',7)
        ylabel('nRPE coefficient (a.u.)', 'FontSize',7)
        hold off;
        title(beta_chan.label{r})
        legend(cplots,names,'FontSize',6,...%'Location','northeast',...
            'NumColumns',1, 'Position',[0.1425,0.7,0.1,0.1])
        %legend('boxoff')
        
        %sgtitle(beta_chan.label{r})
        
    end
    plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa_trajectories.pdf'];
    print(plot_fname,cf,'-dpdf','-fillpage','-r300')
    close all
end