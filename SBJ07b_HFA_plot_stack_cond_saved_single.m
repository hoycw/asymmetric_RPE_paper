function SBJ07b_HFA_plot_stack_cond_saved_single(SBJ, elec_lab, conditions, an_id, actv_win,...
                                        plt_id, save_fig, fig_vis, fig_ftype)
% Plots single trial stack for both stimulus- and response-locked HFA computed in SBJ08a_HFA_actv
%   sorts by condition, then by RT; scatter for RTs in stim-locked
% clear all; %close all;

if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load Results
% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');

% Load data
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
load(hfa_fname,'hfa');

% Load ROI and GM/WM info
% einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' proc_id '.mat'];
% load(einfo_filename);

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim;
cfg_trim.channel = elec_lab;
hfa = ft_selectdata(cfg_trim,hfa);

% Compile cond_type, hit, RT, trl_n
[cond_lab, cond_colors, ~, cond_mrkrs] = fn_condition_label_styles(conditions);
trl_info.cond_n = fn_condition_index(conditions,trl_info);
cond_mat = horzcat(trl_info.cond_n,trl_info.rt,[1:numel(trl_info.trl_n)]');
cond_mat = sortrows(cond_mat,[1 2]);
cond_edges = find(diff(cond_mat(:,1)));
% RT_mean = mean(round(1000*trl_info.rt)); % converts sec to ms
% % Add in the baseline offset to plot correctly
% RT_mean = RT_mean-plt_vars.plt_lim(1)*1000;

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/' SBJ '/stack_' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

% Create a figure for each channel

for ch_ix = 1:numel(hfa.label)
    % Plot parameters
    fig_name = [SBJ '_' conditions '_stack_' hfa.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Get color limits
    clims = NaN([1 2]);
    clims(1) = prctile(reshape(hfa.powspctrm(:,ch_ix,:,:),...
        [1 size(cond_mat,1)*numel(hfa.time)]),plt_vars.clim_perc(1));
    clims(2) = prctile(reshape(hfa.powspctrm(:,ch_ix,:,:),...
        [1 size(cond_mat,1)*numel(hfa.time)]),plt_vars.clim_perc(2));
    clims = [min(clims(1)) max(clims(2))];
    
    %         subplot(1,numel(event_lab),sr_ix);
    hold on;
    ax     = gca;
    
    % Plot Single Trials Per Condition
    imagesc(hfa.time,1:size(cond_mat,1),squeeze(hfa.powspctrm(cond_mat(:,3),ch_ix,:,:)));
    set(gca,'YDir','normal');
    scat = gobjects(size(cond_lab));
    if any(strcmp(plt_vars.evnt_lab,'R'))
        for cond_ix = 1:numel(cond_lab)
            idx = cond_mat(:,1)==cond_ix;
            scat(cond_ix) = scatter(cond_mat(idx,2),find(idx),...
                plt_vars.evnt_mrkr_sz,'MarkerFaceColor',cond_colors{cond_ix},'MarkerEdgeColor','k',...
                'Marker',cond_mrkrs{cond_ix});
        end
    end
    ylim([1 size(cond_mat,1)]);
    
    % Plot events: stim, target, feedback onset, feedback offset
%     plt_vars.evnt_lab    = {'stim','target','fb on','fb off'};
%     plt_vars.evnt_styles = {'-', '--', '-', '-'};
    if numel(plt_vars.evnt_lab)==1
        evnt_times = -plt_vars.plt_lim(1);
    else
        evnt_times = zeros(size(plt_vars.evnt_lab));
        for e = 1:numel(plt_vars.evnt_lab)
            switch plt_vars.evnt_lab{e}
                case 'S'
                    evnt_times(e) = 0;
                case 'R'
                    evnt_times(e) = trl_info.prdm.target;
                case {'Fon','F'}
                    evnt_times(e) = trl_info.prdm.target+trl_info.prdm.fb_delay;
                case 'Foff'
                    evnt_times(e) = trl_info.prdm.trl_len;
                otherwise
                    error('unknown evnt_lab');
            end
        end
    end
    evnt_lines = gobjects(size(evnt_times));
    for evnt_ix = 1:numel(evnt_times)
        evnt_lines(evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle',plt_vars.evnt_styles{evnt_ix});
    end
    
    % Plot condition dividers
    for div_ix = 1:numel(cond_edges)
        line(plt_vars.plt_lim, [cond_edges(div_ix) cond_edges(div_ix)],...
            'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle','--');
    end
    
    % Plotting parameters
    ax = gca;
%     ax.legend        = plt_vars.legend;
    ax.Title.String  = strcat(hfa.label{ch_ix});%, ' (', einfo(ch_ix,2), '): stim trials');
    ax.XLim          = plt_vars.plt_lim;
    ax.XTick         = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Trials';
    cbar = colorbar;
    caxis(clims);
    if plt_vars.legend
        if any(strcmp(plt_vars.evnt_lab,'R'))
            legend(scat,cond_lab{:}, 'Location',plt_vars.legend_loc);
        else
            legend(evnt_lines, plt_vars.evnt_lab, 'Location',plt_vars.legend_loc);
        end
%         legend(event_lines,event_lab,'Location',plt_vars.legend_loc);
    end
    set(gca,'FontSize',16);
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        % Ensure vector graphics if saving
        if any(strcmp(fig_ftype,{'svg','eps'}))
            set(gcf, 'Renderer', 'painters');
        end
        saveas(gcf,fig_fname);
    end
end

end
