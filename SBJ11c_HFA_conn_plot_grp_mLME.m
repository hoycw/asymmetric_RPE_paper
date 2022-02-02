function SBJ11c_HFA_conn_plot_grp_mLME(proc_id, an_id, model_id, conn_id, swap_Xcorr)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

swap_rois = 0;
if strcmp(conn_id,'Xcorr') & swap_Xcorr == 1
    swap_rois = 1;
end
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
[reg_lab, reg_names, colors, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/conn/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

stats_fname = [stats_dir proc_id '_' model_id '_' an_id '_' conn_id '.mat'];
load(stats_fname)

%% plot
ncols = numel(conn_stats.coefs);
if swap_rois == 1
    time_error= [fliplr(conn_stats.time),conn_stats.time];
else
    time_error=[conn_stats.time,fliplr(conn_stats.time)];
end
cf = figure('units','normalized','outerposition',[0 0 1 1],...
            'PaperOrientation','Landscape');
ylims0 = NaN(ncols,2);
for pl = 1:ncols
    cupper = conn_stats.upper{pl};
    clower = conn_stats.lower{pl};
    ylims0(pl,:) = [min(clower(2:end,:),[],'all'), max(cupper(2:end,:),[],'all')];
    ylims0(pl,1) = min([ylims0(pl,1) - 0.2*abs(ylims0(pl,1)), 0]);
    ylims0(pl,2) = max([ylims0(pl,2) + 0.2*abs(ylims0(pl,2)), 0]);
end
ylims = [min(ylims0(:,1)),max(ylims0(:,2))];

for pl = 1:ncols
    ccoef = conn_stats.coefs{pl};
    subplot(1,ncols,pl)
    yline(0); hold on; xline(0); hold on;
    lgnds = NaN(numel(conn_stats.feature),1);
    for rg = 1:numel(conn_stats.feature)
        ccoefrg = ccoef(rg + 1, :);
        coef_error=[conn_stats.lower{pl}(rg + 1,:),fliplr(conn_stats.upper{pl}(rg + 1,:))];
        qmask = double(conn_stats.qvals{pl}(rg + 1,:) < .05);
        qmask(qmask == 0) = Inf;
        if swap_rois == 1
            ccoefrg = fliplr(ccoefrg);
            coef_error = fliplr(coef_error);
            qmask = fliplr(qmask);
        end
        fill(time_error, coef_error,colors{rg},'LineStyle','none'); alpha(0.05); hold on;
        lgnds(rg) = plot(conn_stats.time,ccoefrg, 'Color',colors{rg}); hold on;
        plot(conn_stats.time,ccoefrg.* qmask, 'Color',colors{rg},'LineWidth',3);
        hold on
    end
    ylim(ylims)
    if swap_rois == 1
        title(sprintf('%s to %s',conn_stats.pair_labels{pl}{2},...
        conn_stats.pair_labels{pl}{1}))
    else
    title(sprintf('%s to %s',conn_stats.pair_labels{pl}{1},...
          conn_stats.pair_labels{pl}{2}))
    end
    xlabel('time lag (s)')
    ylabel('model coefficient (a.u.)')
    legend(lgnds,reg_names)  
end
plot_fname = [fig_dir proc_id '_' model_id '_' an_id '_' conn_id '.pdf'];
print(plot_fname,cf,'-dpdf','-fillpage')
close all