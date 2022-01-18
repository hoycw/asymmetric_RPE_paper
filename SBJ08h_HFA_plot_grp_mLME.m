function SBJ08h_HFA_plot_grp_mLME(proc_id, an_id, model_id)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
[reg_lab, reg_names, colors, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa.mat'];
load(stats_fname)

%% plot

ncols = numel(beta.coefs);
time_error=[beta.time,fliplr(beta.time)];
cf = figure('units','normalized','outerposition',[0 0 1 1],...
            'PaperOrientation','Landscape');
ylims0 = NaN(ncols,2);
for r = 1:ncols
    ccoef = beta.coefs{r};
    ylims0(r,:) = [min(ccoef(2:end,:),[],'all'), max(ccoef(2:end,:),[],'all')];
    ylims0(r,1) = min([ylims0(r,1) - 0.5*abs(ylims0(r,1)), 0]);
    ylims0(r,2) = max([ylims0(r,2) + 0.5*abs(ylims0(r,2)), 0]);
end
ylims = [min(ylims0(:,1)),max(ylims0(:,2))];

for r = 1:ncols
    ccoef = beta.coefs{r};
    subplot(1,ncols,r)
    yline(0); hold on; xline(0); hold on;
    lgnds = NaN(numel(reg_lab),1);
    for rg = 1:numel(reg_lab)
        coef_error=[beta.lower{r}(rg + 1,:),fliplr(beta.upper{r}(rg + 1,:))];
        qmask = double(beta.qvals{r}(rg + 1,:) < .05);
        qmask(qmask == 0) = Inf;
        fill(time_error, coef_error,colors{rg},'LineStyle','none'); alpha(0.05); hold on;
        lgnds(rg) = plot(beta.time,ccoef(rg + 1,:), 'Color',colors{rg}); hold on;
        plot(beta.time,ccoef(rg + 1,:).* qmask, 'Color',colors{rg},'LineWidth',3);
        hold on
    end
    ylim(ylims)
    title(sprintf('%s',beta.label{r}))
    xlabel('time lag (s)')
    ylabel('model coefficient (a.u.)')
    legend(lgnds,reg_names)  
end
plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_hfa.pdf'];
print(plot_fname,cf,'-dpdf','-fillpage')
close all