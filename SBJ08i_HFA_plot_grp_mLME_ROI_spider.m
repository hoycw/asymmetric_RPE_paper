function SBJ08i_HFA_plot_grp_mLME_ROI_spider(SBJ_id,proc_id,an_id,model_id,stat_id,cat_id,...
                                          roi_id,save_fig,varargin)
%% Plot spider plot of proportion of electrodes per puns category for each ROI
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   model_id [str] - ID of the model used in GLM
%   stat_id [str] - ID of the statistical parameters to extract from HFA
%   roi_id [str] - ID of set of ROIs to group and color electrodes
%   plot_scat [0/1] - binary flag to plot scatter points for individual SBJs
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       atlas_id [str] - ID of the atlas for elec ROI assignments
%           default: 'Dx'
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   Bar chart with proprotion of ROI showing each effect

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir();
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);

%% Handle Variable Inputs & Defaults
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

%% Prep variables
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
if ~strcmp(cat_id,'puns'); error('not ready for non-puns categories'); end

SBJs = fn_load_SBJ_list(SBJ_id);

% Get condition info
[reg_lab, reg_names, reg_colors, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
[cat_lab, cat_names, cat_colors, ~, ~] = fn_puns_category_label_styles(cat_id);

% Load all ROI info
[roi_list, roi_colors, roi_field] = fn_roi_label_styles(roi_id);

% Set up electrode counts
reg_cnt  = zeros([numel(SBJs) numel(roi_list) numel(reg_lab)]);
elec_cnt = zeros([numel(SBJs) numel(roi_list)]);
out_cnt  = zeros(size(SBJs));
miss_cnt = zeros(size(SBJs));

%% Load stats results:
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/proportions/' model_id '/' stat_id '/' an_id '/' atlas_id '_' roi_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [root_dir 'PRJ_Error/data/GRP/stats/' model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname,'beta_chan');

%% Load significant channel labels
roi_cat_elecs = cell(numel(roi_list),numel(cat_lab));
for roi_ix = 1:numel(roi_list)
    beta_roi_ix = strcmp(beta_chan.label,roi_list{roi_ix});
    for cat_ix = 1:length(cat_lab)
        beta_cat_ix = strcmp(beta_chan.chancat_label,cat_lab{cat_ix});
        elec_ix = beta_chan.chancat_ix{beta_roi_ix}{beta_cat_ix,1};
        roi_cat_elecs{roi_ix,cat_ix} = beta_chan.chan_label{beta_roi_ix}(elec_ix);
    end
end

%% Load Electrode Counts and Compute Within-Subject Proportions
roi_cat_prop = zeros(numel(roi_list), numel(cat_lab), numel(SBJs));
roi_sbj = true(numel(roi_list), numel(SBJs));
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    
    % Aggregate results per ROI
    for roi_ix = 1:numel(roi_list)
        roi_ch_idx = strcmp(elec.(roi_field),roi_list{roi_ix});
        if any(roi_ch_idx)
            for cat_ix = 1:numel(cat_lab)
                sbj_elecs = split(roi_cat_elecs{roi_ix,cat_ix}(contains(roi_cat_elecs{roi_ix,cat_ix},SBJ)),' ');
                if ~isempty(sbj_elecs)
                    % flip in case single electrode comes out as column vector
                    if size(sbj_elecs,2)==1; sbj_elecs = sbj_elecs'; end
                    sbj_elecs = sbj_elecs(:,2); % toss SBJ column
                    % Double check all significant electrodes are in elec
                    for e_ix = 1:numel(sbj_elecs)
                        if ~any(strcmp(elec.label(roi_ch_idx),sbj_elecs{e_ix}))
                            error([SBJ ' elec ' sbj_elecs{e_ix} ' not found in '...
                                roi_list{roi_ix} ' in elec file!']);
                        end
                    end
                    % Compute proportion
                    roi_cat_prop(roi_ix,cat_ix,sbj_ix) = numel(sbj_elecs)/sum(roi_ch_idx);
                end
            end
        else
            roi_sbj(roi_ix,sbj_ix) = false;
            roi_cat_prop(roi_ix,:,sbj_ix) = nan;
        end
    end
    clear SBJ SBJ_vars elec
end

% Compute mean and SEM
roi_cat_mean = nanmean(roi_cat_prop,3);
for roi_ix = 1:numel(roi_list)
    roi_cat_sem(roi_ix,:)  = nanstd(roi_cat_prop(roi_ix,:,:),[],3)./sqrt(sum(roi_sbj(roi_ix,:)));
end

%% Plot Proportion of Electrodes by Category and ROI
% Create and format the plot
fig_name = [SBJ_id '_HFA_spider_' model_id '_' stat_id '_' cat_id '_' atlas_id '_' roi_id];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

cat_order = {'uRPE','pRPE','sRPE','nRPE'};
for cat_ix = 1:numel(cat_lab)
    cat_idx(cat_ix) = find(strcmp(cat_lab,cat_order{cat_ix}));
end
max_lim = max(roi_cat_mean(:)) - mod(max(roi_cat_mean(:)),0.05) + 0.05;
interval = round(max_lim/0.05);
spider = spider_plot(roi_cat_mean(:,cat_idx),'AxesLabels',cat_names(cat_idx),...
    'Color',vertcat(roi_colors{:}),'LineWidth',4,'MarkerSize',75,...
    'AxesLimits',[zeros(1,numel(cat_lab)); ones(1,numel(cat_lab))*max_lim],...
    'AxesInterval',interval,'AxesPrecision',2,'AxesFontSize',14,'LabelFontSize',18);

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Bars
if plot_scat
    scat_suffix = '_SBJscat';
else
    scat_suffix = '';
end
fig_name = [SBJ_id '_HFA_errbar_' model_id '_' stat_id '_' cat_id '_' atlas_id '_' roi_id scat_suffix];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% Plot Category bars within each ROI
bars = cell(size(roi_list));
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(cat_lab));
for roi_ix = 1:numel(roi_list)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    bars{roi_ix} = bar(bar_offsets+roi_ix,diag(roi_cat_mean(roi_ix,:)),0.9,'stacked');
    for cat_ix = 1:numel(cat_lab)
        set(bars{roi_ix}(cat_ix),'FaceColor',cat_colors{cat_ix},'EdgeColor','k');
        
        % Plot standard error of the mean across subjects
        line([bar_offsets(cat_ix) bar_offsets(cat_ix)]+roi_ix,...
             [roi_cat_mean(roi_ix,cat_ix)+roi_cat_sem(roi_ix,cat_ix) ...
              roi_cat_mean(roi_ix,cat_ix)-roi_cat_sem(roi_ix,cat_ix)],...
            'Color','k','LineWidth',1.5);
        
        % Overlay scatter for individual SBJ
        if plot_scat
            scat = scatter(scat_offsets(roi_sbj(roi_ix,:))+bar_offsets(cat_ix)+roi_ix,...
                roi_cat_prop(roi_ix,cat_ix,roi_sbj(roi_ix,:)),50,'k*');
        end
    end
end
leg_loc = 'northwest';
if plot_scat
    legend([bars{1},scat],cat_names{:},'Individuals','Location',leg_loc);
else
    legend(bars{1},cat_names{:},'Location',leg_loc);
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim       = [0.5 0.5+numel(roi_list)];
ax.XTick      = 1:numel(roi_list);
ax.XColor     = 'k';
ax.XTickLabel = roi_list;

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLim            = [0 ax.YLim(2)];%[0 0.6];%
% ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Proportion of Electrodes per RPE Category and ROI';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
set(gca,'FontSize',16);

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
