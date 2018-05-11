function SBJ08b_HFA_plot_stack_cond_saved(SBJ,conditions,pipeline_id, an_id,actv_win,...
                                        plt_id,save_fig,fig_vis,fig_filetype)
% Plots single trial stack for both stimulus- and response-locked HFA computed in SBJ08a_HFA_actv
%   sorts by condition, then by RT; scatter for RTs in stim-locked
% clear all; %close all;

if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Error/scripts/');
addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');

% Load data
hfa_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id,'.mat');
load(hfa_filename,'hfa');

% Load ROI and GM/WM info
einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
load(einfo_filename);
% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa.time)-1)/(hfa.time(end)-hfa.time(1));

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim;
hfa = ft_selectdata(cfg_trim,hfa);

% Compile cond_type, hit, RT, trl_n
[cond_lab, cond_colors, ~, cond_mrkrs] = fn_condition_label_styles(conditions);
trl_info.cond_n = fn_condition_index(conditions,trl_info);
cond_mat = horzcat(trl_info.cond_n,round(1000*trl_info.rt),[1:numel(trl_info.trl_n)]');
cond_mat = sortrows(cond_mat,[1 2]);
cond_edges = find(diff(cond_mat(:,1)));
RT_mean = mean(round(1000*trl_info.rt)); % converts sec to ms
% Add in the baseline offset to plot correctly
RT_mean = RT_mean-plt_vars.plt_lim(1)*1000;

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Error/results/HFA/' SBJ '/stack_' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa.label)
    % Plot parameters
    fig_name = [SBJ '_' conditions '_stack_' hfa.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
    
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
    imagesc(squeeze(hfa.powspctrm(cond_mat(:,3),ch_ix,:,:)));
    set(gca,'YDir','normal');
    x_tick_lab = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    scat = [];
    for cond_ix = 1:numel(cond_lab)
        idx = cond_mat(:,1)==cond_ix;
        scat(cond_ix) = scatter(cond_mat(idx,2)-plt_vars.plt_lim(1)*sample_rate,find(idx),...
            'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k',...
            'Marker',cond_mrkrs{cond_ix});
    end
    ylim([1 size(cond_mat,1)]);
    
    % Plot events: stim, target, feedback onset, feedback offset
    event_lab    = {'stim','target','fb on','fb off'};
    event_styles = {'-', '--', '-', '-'};
    event_times = [-plt_vars.plt_lim(1)*sample_rate...
        (trl_info.timing.target-plt_vars.plt_lim(1))*sample_rate...
        (trl_info.timing.target+trl_info.timing.fb_delay-plt_vars.plt_lim(1))*sample_rate...
        (trl_info.timing.trl_len-plt_vars.plt_lim(1))*sample_rate];
    event_lines = [];
    for event_ix = 1:numel(event_times)
        event_lines(event_ix) = line([event_times(event_ix) event_times(event_ix)],ylim,...
            'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle',event_styles{event_ix});
    end
    
    % Plot condition dividers
    for div_ix = 1:numel(cond_edges)
        line([0,numel(hfa.time)], [cond_edges(div_ix) cond_edges(div_ix)],...
            'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle','--');
    end
    
    % Plotting parameters
    ax = gca;
%     ax.legend        = plt_vars.legend;
    ax.Title.String  = strcat(hfa.label{ch_ix}, ' (', einfo(ch_ix,2), '): stim trials');
    ax.XLim          = [0,numel(hfa.time)];
    ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:numel(hfa.time);
    ax.XTickLabel    = x_tick_lab;
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Trials';
    cbar = colorbar;
    caxis(clims);
    if plt_vars.legend
        legend(scat,cond_lab{:}, 'Location',plt_vars.legend_loc);
%         legend(event_lines,event_lab,'Location',plt_vars.legend_loc);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
