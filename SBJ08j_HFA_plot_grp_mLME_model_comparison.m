function SBJ08j_HFA_plot_grp_mLME_model_comparison(model_ids, model_names, an_id, stat_id)
%% set directories
if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/'; ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
%% Load models
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];
crit = [];
crit.AIC = {}; crit.MSE = {}; crit.r2 = {}; crit.BIC = {};
for mp = 1:numel(model_ids)
    beta = [];
    cmodel= model_ids{mp};
    display(cmodel)
    load([stats_dir cmodel '_' stat_id '_' an_id '_hfa.mat']);
    display(beta.feature)
    for r = 1:numel(beta.AIC)
        crit.AIC{r}(mp,:) = beta.AIC{1,r};
        crit.MSE{r}(mp,:) = beta.MSE{1,r};
        crit.BIC{r}(mp,:) = beta.BIC{1,r};
        crit.r2{r}(mp,:) = beta.r2{1,r};
    end
    if mp == 1
        time = beta.time;
        label = beta.label;
    end
end
%% Make a plot
% cf1 = figure('units','normalized','outerposition',[0 0 1 1],...
%     'PaperOrientation','Landscape');
% criteria = fieldnames(crit);
% for c = 1:numel(criteria)
%     ccrit = criteria{c};
%     for r = 1:numel(crit.(ccrit))
%         subplot(ceil(numel(crit.(ccrit)) / 2), 2, r)
%         plot(time, crit.(ccrit){r}(2,:) - crit.(ccrit){r}(1,:)); hold on;
%     end
% end
%
%% Make a plot
%model_names = cellfun(@(x) x(1:end-6), model_ids,'UniformOutput',false);
criteria = fieldnames(crit);
markers = 'dhs*';
lty = {'-.','--','-',':'};
for c = 1:numel(criteria)
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperOrientation','Landscape');
    ccrit = criteria{c};
    for r = 1:numel(crit.(ccrit))
        subplot(ceil(numel(crit.(ccrit)) / 2), 2, r)
        for m = 1:size(crit.(ccrit){r},1)
            cformat = lty{m};% markers(m)];
            cp = plot(time, crit.(ccrit){r}(m,:), cformat, 'LineWidth',1.5); hold on;
            cp.Color(4) = 0.8;
        end
        title(label{r})
        xlabel('time')
        ylabel(ccrit)
    end
    sgtitle(['Model comparison (' ccrit ')'])
    legend(model_names,'Location','southeast')
    plot_fname = [fig_dir stat_id '_' an_id '_hfa_model_comparison_' ccrit '.pdf'];
    print(plot_fname,cf,'-dpdf','-fillpage')
end
close all
