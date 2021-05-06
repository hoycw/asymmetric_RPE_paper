function SBJ08d_HFA_plot_grp_GLM_ts_gROIcomb(SBJ_id,proc_id,an_id,model_id,stat_id,...
                                    roi_id,plt_id,save_fig,varargin)
% Plots mGLM beta time series by ROI groupings
%   gROI version: one plot; subplots per gROI (lines still colored by ROI)
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of statistical analysis
%   an_id [str] - analysis ID for S-locked preprocessing, filtering, etc.
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   plt_id [str] - ID of the plotting variables
%   save_fig [0/1] - save this figure?
%   fig_vis [str] - visible = 'on'/'off'
%   fig_ftype [str] - file extension for figure saving

%   gm_thresh [float] - threshold of GM % to include electrode (likely = 0)
%   z_thresh [float] - threshold of HFA z score to include electrode
%   plot_nsig [0/1] - plot electrodes with no significant epochs?

%% Handle Variable Inputs & Defaults
% add %gm_thresh,z_thresh,plot_nsig,...
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        elseif strcmp(varargin{v},'atlas_id') && ischar(varargin{v+1})
            atlas_id = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');   fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('atlas_id','var');  atlas_id = 'Dx'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);
if ~strcmp(plt.sig_type,'bold')
    error('Only bold sig_type allowed');
end

SBJs = fn_load_SBJ_list(SBJ_id);

% Get condition info
[reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'gROI','mgROI'}))
    roi_field = 'gROI';
else
    error('not ready for non gROI');
end


%% Load Data
% Load example to get timing info
load([root_dir 'PRJ_Error/data/' SBJs{1} '/03_events/' SBJs{1} '_bhv_' proc_id '_final.mat'],'bhv');
load([root_dir 'PRJ_Error/data/' SBJs{1} '/07_stats/' SBJs{1} '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat']);
time_vec = beta.time;

% Get event times for plotting
% [evnt_times] = fn_get_evnt_times(an.evnt_lab,plt.evnt_lab,bhv);

sig_ts  = zeros([numel(reg_lab) numel(roi_list) numel(time_vec)]);
sig_ch  = cell(size(SBJs));
roi_ch  = cell(size(SBJs));
onsets  = cell([numel(reg_lab) numel(roi_list)]);
roi_cnt = zeros(size(roi_list));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([SBJ_vars.dirs.recon SBJs{s} '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    load([SBJ_vars.dirs.stats SBJs{s} '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat']);
    
    sig_ch{s} = false([numel(beta.label) numel(reg_lab)]);
    roi_ch{s} = zeros(size(beta.label));
    for ch_ix = 1:numel(beta.label)
        elec_ix = strcmp(beta.label{ch_ix},elec.label);
        if any(elec_ix)
            roi_ix = find(strcmp(roi_list,elec.(roi_field){elec_ix}));
            if ~isempty(roi_ix)
                roi_ch{s}(ch_ix) = roi_ix;
                roi_cnt(roi_ix) = roi_cnt(roi_ix) + 1;
                % Consolidate to binary sig/non-sig
                for reg_ix = 1:numel(reg_lab)
                    beta_reg_ix = strcmp(beta.feature,reg_lab{reg_ix});
                    if any(beta.qval(beta_reg_ix,ch_ix,:)<=0.05,3)
                        sig_ch{s}(ch_ix,reg_ix) = true;
                        sig_ts(reg_ix,roi_ix,:) = squeeze(sig_ts(reg_ix,roi_ix,:))...
                                                + squeeze(beta.qval(beta_reg_ix,ch_ix,:)<=st.alpha);
                        onsets{reg_ix,roi_ix} = [onsets{reg_ix,roi_ix}...
                            min(time_vec(beta.qval(beta_reg_ix,ch_ix,:)<=st.alpha))];
                    end
                end
            end
        else
            fprintf(2,'\t%s: %s in beta missing from elec!\n',SBJ,beta.label{ch_ix});
        end
    end
    
    clear beta elec SBJ SBJ_vars SBJ_vars_cmd
end

% % HACK!!! to mesh with PRJ_Stroop elec files
% cfgs = [];
% if strcmp(SBJ,'CP24')
%     cfgs.channel = {'all','-RTO4'};
% elseif strcmp(SBJ,'IR57')
%     cfgs.channel = {'all','-LAM4-5','-LAM5-6','-RIN8-9','-RIN9-10',...
%         '-RSM1-2','-RSM2-3','-RSM3-4','-RTI9-10'};
% elseif strcmp(SBJ,'IR68')
%     cfgs.channel = {'all','-LPC5-6','-LPC6-7','-LPC7-8'};
% end
% stat = ft_selectdata(cfgs,stat);
% beta   = ft_selectdata(cfgs,beta);

% % Sort elecs by stat labels
% if ~strcmp(SBJ,'IR66')  % HACK!!!
%     cfgs = []; cfgs.channel = stat.label;
%     elec = fn_select_elec(cfgs,elec);
%     elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
% else
%     elec.roi = elec.(roi_field);
% end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP_GLM_ts/' model_id '/' stat_id '/' an_id '/gROIcomb/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

ylims = [-0.5 max(sig_ts(:))+1];

% Create a figure for each condition (all gROIs in one subplot)
fig_name = [SBJ_id '_GLM_ts_' model_id '_' stat_id '_' atlas_id '_' roi_id];
f = figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis); %twice as wide for the double plot
for reg_ix = 1:numel(reg_lab)
    ax = subplot(numel(reg_lab),1,reg_ix);%subplot_ix);
    hold on;
    roi_lines = gobjects(size(roi_list));
    for roi_ix = 1:numel(roi_list)
        roi_lines(roi_ix) = plot(time_vec,squeeze(sig_ts(reg_ix,roi_ix,:)),...
            'Color',roi_colors{roi_ix},'LineWidth',2);%,'LineStyle',grp_style{cond_ix});
    end
    
%     evnt_lines = gobjects(size(evnt_times));
%     for evnt_ix = 1:numel(evnt_times)
%         evnt_lines(evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
%             'LineWidth',plt.evnt_width,'Color','k','LineStyle',plt.evnt_styles{evnt_ix});
%     end
    
    % Plotting parameters
    ax.Title.String  = reg_names{reg_ix};
    ax.Box           = 'off';
    ax.YLim          = ylims;
%     ax.YTick         = yticks;
    ax.YLabel.String = '# Sig Elecs';
%     ax.XLim          = [eval(['plt.plt_lim_' evnt_labs]);%[min(beta.win_lim_s(:,1)) max(beta.win_lim_s(:,2))];
    %                 ax.XTick         = 0:plt.x_step_sz*srate:size(beta.time,2);
%     ax.XTick         = x_tick_lab;
    ax.XLabel.String = 'Time (s)';
    
    % Legend
    roi_legend = cell(size(roi_list));
    for roi_ix = 1:numel(roi_list)
        roi_legend{roi_ix} = [roi_list{roi_ix} ' (n = ' num2str(roi_cnt(roi_ix)) ')'];
    end
    legend(roi_lines,roi_legend,'Location',plt.legend_loc);
    
    set(gca,'FontSize',16);
end

%% Print onset latency stats
for reg_ix = 1:numel(reg_lab)
    for roi_ix = 1:numel(roi_list)
        if ~isempty(onsets{reg_ix,roi_ix})
            fprintf('%s %s: n = %d / %d (%.03f); %.03f +/- %.04f ms (min = %.03f, max = %.03f)\n', reg_lab{reg_ix}, roi_list{roi_ix},...
                numel(onsets{reg_ix,roi_ix}), roi_cnt(roi_ix), numel(onsets{reg_ix,roi_ix})/roi_cnt(roi_ix), ...
                mean(onsets{reg_ix,roi_ix}), std(onsets{reg_ix,roi_ix}),...
                min(onsets{reg_ix,roi_ix}), max(onsets{reg_ix,roi_ix}));
        else
            fprintf(2,'%s %s: n = 0\n', reg_lab{reg_ix}, roi_list{roi_ix});
        end
    end
end

%% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
