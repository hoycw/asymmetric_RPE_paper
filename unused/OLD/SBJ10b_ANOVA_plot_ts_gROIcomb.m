function SBJ10b_ANOVA_plot_ts_gROIcomb(SBJ,proc_id,stat_id,an_id,...
                                    atlas_id,roi_id,gm_thresh,z_thresh,plot_nsig,...
                                    plt_id,save_fig,fig_vis,fig_ftype)
% Plots smANOVA w2 time series by ROI groupings
%   gROI version: one plot; subplots per gROI (lines still colored by ROI)
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of statistical analysis
%   an_id [str] - analysis ID for S-locked preprocessing, filtering, etc.
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   gm_thresh [float] - threshold of GM % to include electrode (likely = 0)
%   z_thresh [float] - threshold of HFA z score to include electrode
%   plot_nsig [0/1] - plot electrodes with no significant epochs?
%   plt_id [str] - ID of the plotting variables
%   save_fig [0/1] - save this figure?
%   fig_vis [str] - visible = 'on'/'off'
%   fig_ftype [str] - file extension for figure saving
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
if ~strcmp(plt_vars.sig_type,'bold')
    error('Only bold sig_type allowed');
end

% Get condition info
[grp_lab, ~, ~] = fn_group_label_styles(model_lab);

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'gROI','mgROI'}))
    roi_field = 'gROI';
else
    error('not ready for non gROI');
end

%% Load Data
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');
load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);

load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',stat_id,'_',an_id,'.mat'));
% Convert to Percentages
w2.trial = w2.trial.*100;

% HACK!!! to mesh with PRJ_Stroop elec files
cfgs = [];
if strcmp(SBJ,'CP24')
    cfgs.channel = {'all','-RTO4'};
elseif strcmp(SBJ,'IR57')
    cfgs.channel = {'all','-LAM4-5','-LAM5-6','-RIN8-9','-RIN9-10',...
        '-RSM1-2','-RSM2-3','-RSM3-4','-RTI9-10'};
elseif strcmp(SBJ,'IR68')
    cfgs.channel = {'all','-LPC5-6','-LPC6-7','-LPC7-8'};
end
stat = ft_selectdata(cfgs,stat);
w2   = ft_selectdata(cfgs,w2);

% Sort elecs by stat labels
if ~strcmp(SBJ,'IR66')  % HACK!!!
    cfgs = []; cfgs.channel = stat.label;
    elec = fn_select_elec(cfgs,elec);
    elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
else
    elec.roi = elec.(roi_field);
end

%% Select data
% FDR correct pvalues for ANOVA
qvals = NaN(size(w2.pval));
for e_ix = 1:numel(w2.label)
    [~, ~, ~, qvals(:,e_ix,:)] = fdr_bh(squeeze(w2.pval(:,e_ix,:)));%,0.05,'pdep','yes');
end

% Find significant channels
sig_idx = false([numel(grp_lab) numel(w2.label)]);
for grp_ix = 1:numel(grp_lab)
    for e_ix = 1:numel(w2.label)
        if any(squeeze(qvals(grp_ix,e_ix,:))<=0.05)
            sig_idx(grp_ix,e_ix) = true;
        end
    end
end
sig_elecs = w2.label(squeeze(any(sig_idx,1)));

% Exclude elecs not in atlas ROIs
roi_elecs = fn_select_elec_lab_match(elec, 'b', atlas_id, roi_id);

% GM Thresh
% Z Thresh

% Select significant, ROI-matched elecs
plot_elecs = intersect(sig_elecs, roi_elecs);
sig_elec_ix = zeros(size(plot_elecs));
for e_ix = 1:numel(plot_elecs)
    sig_elec_ix(e_ix) = find(strcmp(w2.label,plot_elecs(e_ix)));
end

cfgs = []; cfgs.channel = plot_elecs;
elec = fn_select_elec(cfgs, elec);
%stat = ft_selectdata(cfgs,stat);
w2   = ft_selectdata(cfgs,w2);
sig_idx = sig_idx(:,sig_elec_ix);
qvals = qvals(:,sig_elec_ix,:);

% Add ROI plotting color
elec.color  = cell(size(elec.label));
for e_ix = 1:numel(elec.label)
    elec.color(e_ix) = roi_colors(strcmp(elec.(roi_field){e_ix},roi_list));
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/' SBJ '/ANOVA_ts/' stat_id '/' an_id '/gROIcomb/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

% Find plot limits
max_w2 = max(max(max(w2.trial)));
min_w2 = min(min(min(w2.trial)));
ylim_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims  = [min_w2-ylim_fudge max_w2+ylim_fudge];
yticks = 0:1:ylims(2);

% Create a figure for each condition (all gROIs in one subplot)
fig_name = [SBJ '_ANOVA_ts_' model_lab '_' atlas_id '_' roi_id];
f = figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis); %twice as wide for the double plot
for grp_ix = 1:numel(grp_lab)
    any_plot = 0;
    ax = subplot(numel(grp_lab),1,grp_ix);%subplot_ix);
    hold on;
    % Subplots for each group
    sig_roi_ix    = find(sig_idx(grp_ix,:));
    grp_roi_list  = unique(elec.(roi_field)(sig_roi_ix));
    roi_sig_count = zeros(size(grp_roi_list));
    roi_lines     = gobjects(size(grp_roi_list));
    for e_ix = sig_roi_ix
        % Plot significant time periods
        any_plot = 1;
        roi_line = plot(w2.time,squeeze(w2.trial(grp_ix,e_ix,:))',...
            'Color',elec.color{e_ix});%,'LineStyle',grp_style{cond_ix});
        
        % Track number of elecs per ROI and legend
        roi_ix = strcmp(grp_roi_list,elec.(roi_field){e_ix});
        roi_sig_count(roi_ix) = roi_sig_count(roi_ix) + 1;
        roi_lines(roi_ix) = roi_line;
        
        % Plot significant epochs in bold
        sig_chunks = fn_find_chunks(squeeze(qvals(grp_ix,e_ix,:))<=0.05);
        sig_chunks(squeeze(qvals(grp_ix,e_ix,sig_chunks(:,1)))>0.05,:) = [];
        for sig_ix = 1:size(sig_chunks,1)
            if sig_chunks(sig_ix,1)==sig_chunks(sig_ix,2)
                % Handle case of single significant window
                scatter(w2.time(sig_chunks(sig_ix,1)),squeeze(w2.trial(grp_ix,e_ix,sig_chunks(sig_ix,1))),...
                    50,elec.color{e_ix},'o','filled');
            else
                line(w2.time(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                    squeeze(w2.trial(grp_ix,e_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                    'Color',[elec.color{e_ix}],...'LineStyle',plt_vars.sig_style,
                    'LineWidth',plt_vars.sig_width);
            end
        end
    end
    
    % Plot events
    if numel(plt_vars.evnt_lab)==1
        evnt_times = -plt_vars.plt_lim(1);
    else
        evnt_times = zeros(size(plt_vars.evnt_lab));
        for e = 1:numel(plt_vars.evnt_lab)
            switch plt_vars.evnt_lab{e}
                case 'S'
                    evnt_times(e) = -plt_vars.plt_lim(1);
                case 'R'
                    evnt_times(e) = (trl_info.prdm.target-plt_vars.plt_lim(1));
                case {'Fon','F'}
                    evnt_times(e) = (trl_info.prdm.target+trl_info.prdm.fb_delay-plt_vars.plt_lim(1));
                case 'Foff'
                    evnt_times(e) = (trl_info.prdm.trl_len-plt_vars.plt_lim(1));
                otherwise
                    error('unknown evnt_lab');
            end
        end
    end
    evnt_lines = gobjects(size(evnt_times));
    for evnt_ix = 1:numel(evnt_times)
        evnt_lines(evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylims,...
            'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle',plt_vars.evnt_styles{evnt_ix});
    end
    
    % Plotting parameters
    ax.Title.String  = [grp_lab{grp_ix} ' (n sig=' num2str(sum(roi_sig_count)) '; '...
        strjoin(grp_roi_list,',') ')'];
    ax.Box           = 'off';
    ax.YLim          = ylims;
%     ax.YTick         = yticks;
    ax.YLabel.String = 'Omega^2';
%     ax.XLim          = [eval(['plt_vars.plt_lim_' evnt_labs]);%[min(w2.win_lim_s(:,1)) max(w2.win_lim_s(:,2))];
    %                 ax.XTick         = 0:plt_vars.x_step_sz*srate:size(w2.time,2);
%     ax.XTick         = x_tick_lab;
    ax.XLabel.String = 'Time (s)';
    
    % Legend
    roi_legend = cell(size(grp_roi_list));
    for roi_ix = 1:numel(grp_roi_list)
        roi_legend{roi_ix} = [roi_list{roi_ix} ' (n sig=' num2str(roi_sig_count(roi_ix)) ')'];
    end
    legend(roi_lines,roi_legend,'Location',plt_vars.legend_loc);
    
    set(gca,'FontSize',16);
end

%% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
