function SBJ08b_HFA_plot_crRT_nANOVA_single(SBJ, elec_lab, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype)
% Plots ANOVA results
% clear all; %close all;
% fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%% Load Results
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

[grp_lab, grp_colors, grp_style] = fn_group_label_styles(st.model_lab);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');
load([SBJ_vars.dirs.proc SBJ '_nANOVA_ROI_' stat_id '_' an_id '.mat']);
% Convert % explained variance to 0-100 scale
w2.trial = w2.trial*100;

% Load ROI and GM/WM info
% einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' proc_id '.mat'];
% load(einfo_filename);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/' SBJ '/' stat_id '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Find plot limits
max_w2 = max(max(max(w2.trial)));
min_w2 = min(min(min(w2.trial)));
ylim_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims  = [min_w2-ylim_fudge max_w2+ylim_fudge];

% y_sig = zeros([1 numel(grp_lab)+1]);
% y_sig(1) = mean([min_w2,max_w2]);
% for grp_ix = 2:numel(grp_lab)+2
%     y_sig(grp_ix) = y_sig(grp_ix-1)+ylim_fudge;
% end

% Create a figure for each channel
if numel(elec_lab)>1; error('make this work for more than one'); end
ch_ix = find(strcmp(w2.label,elec_lab));

fig_name = [SBJ '_' stat_id '_' w2.label{ch_ix}];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 0.8],'Visible',fig_vis);

ax = gca; hold on;
main_lines = [];
% Plot var_exp
for grp_ix = 1:numel(grp_lab)
    main_lines(grp_ix) = plot(w2.time,squeeze(w2.trial(grp_ix,ch_ix,:)),...
            'Color',grp_colors{grp_ix},'LineStyle',grp_style{grp_ix});
end

% Plot events: stim, target, feedback onset, feedback offset
%     plt_vars.evnt_lab    = {'stim','target','fb on','fb off'};
%     plt_vars.evnt_styles = {'-', '--', '-', '-'};
% if numel(plt_vars.evnt_lab)==1
%     evnt_times = -plt_vars.plt_lim(1)*sample_rate;
% else
%     evnt_times = zeros(size(plt_vars.evnt_lab));
%     for e = 1:numel(plt_vars.evnt_lab)
%         switch plt_vars.evnt_lab{e}
%             case 'S'
%                 evnt_times(e) = -plt_vars.plt_lim(1)*sample_rate;
%             case 'R'
%                 evnt_times(e) = (trl_info.prdm.target-plt_vars.plt_lim(1))*sample_rate;
%             case {'Fon','F'}
%                 evnt_times(e) = (trl_info.prdm.target+trl_info.prdm.fb_delay-plt_vars.plt_lim(1))*sample_rate;
%             case 'Foff'
%                 evnt_times(e) = (trl_info.prdm.trl_len-plt_vars.plt_lim(1))*sample_rate;
%             otherwise
%                 error('unknown evnt_lab');
%         end
%     end
% end
% evnt_lines = [];
% for evnt_ix = 1:numel(evnt_times)
%     evnt_lines(evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylims,...
%         'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle',plt_vars.evnt_styles{evnt_ix});
% end
%         % Plot Response Marker
%         event_line = line([find(stat.time==0) find(stat.time==0)],ylims1,...
%             'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%             'LineStyle',plt_vars.evnt_style);
%         main_lines = [main_lines event_line];

% Plot significant time periods
for grp_ix = 1:numel(grp_lab)
    if any(squeeze(w2.qval(grp_ix,ch_ix,:))<=st.alpha)
        % Find significant periods
        if strcmp(plt_vars.sig_type,'bold')
            sig_chunks = fn_find_chunks(squeeze(w2.qval(grp_ix,ch_ix,:))<=st.alpha);
            sig_chunks(squeeze(w2.qval(grp_ix,ch_ix,sig_chunks(:,1)))>st.alpha,:) = [];
            for sig_ix = 1:size(sig_chunks,1)
                if diff(sig_chunks(sig_ix,:))==0
                    scatter(w2.time(sig_chunks(sig_ix,1)),squeeze(w2.trial(grp_ix,ch_ix,sig_chunks(sig_ix,1))),...
                        plt_vars.sig_scat_size,grp_colors{grp_ix},plt_vars.sig_scat_mrkr,'filled');
                else
                    line(w2.time(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                        squeeze(w2.trial(grp_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                        'Color',grp_colors{grp_ix},'LineStyle',plt_vars.sig_style,...
                        'LineWidth',plt_vars.sig_width);
                end
            end
            %                 elseif strcmp(plt_vars.sig_type,'scatter')
            %                     sig_times = win_center(squeeze(qvals(grp_ix,ch_ix,:))<0.05);
            %                     scatter(sig_times,repmat(y_sig(grp_ix),size(sig_times)),...
            %                         plt_vars.sig_scat_size,grp_colors{grp_ix},plt_vars.sig_scat_mrkr);
        elseif strcmp(plt_vars.sig_type,'patch')
            sig_chunks = fn_find_chunks(squeeze(qvals(grp_ix,ch_ix,:))<0.05);
            sig_chunks(squeeze(qvals(grp_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
            error('sig_type = patch needs sig_times variable!');
            % Plot Significance Shading
            fprintf('%s %s -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                stat.label{ch_ix},grp_lab{grp_ix},size(sig_chunks,1));
            for sig_ix = 1:size(sig_chunks,1)
                sig_times = win_lim(sig_chunks(sig_ix,:),:);
                patch([sig_times(1,1) sig_times(1,1) sig_times(2,2) sig_times(2,2)], ...
                    [ylims(1) ylims(2) ylims(2) ylims(1)],...
                    grp_colors{grp_ix},'FaceAlpha',plt_vars.sig_alpha);
            end
        end
    else
        fprintf('%s %s -- NO SIGNIFICANT CLUSTERS FOUND...\n',...
            w2.label{ch_ix},grp_lab{grp_ix});
    end
end

% Plotting parameters
ax.Title.String  = strcat(w2.label{ch_ix});%, ' (', einfo(ch_ix,2), ')');
ax.Box           = 'off';
% ax.YLim          = ylims;
ax.YLabel.String = '% Variance Explained';
ax.XLim          = plt_vars.plt_lim;
ax.XTick         = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
ax.XLabel.String = 'Time (s)';
legend(main_lines,grp_lab,'Location',plt_vars.legend_loc);
set(ax,'FontSize',16);

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
    %eval(['export_fig ' fig_filename]);
end

end
