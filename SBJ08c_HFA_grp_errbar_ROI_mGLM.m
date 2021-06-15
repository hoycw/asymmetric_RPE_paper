function SBJ08c_HFA_grp_errbar_ROI_mGLM(SBJ_id,proc_id,an_id,model_id,stat_id,roi_id,...
                                          plot_scat,save_fig,varargin)
%% Plot bar graph of proportion of significant electrodes for mass GLM
%   Future: may also plot activation and RT correlations
% INPUTS:
%   SBJ [str] - ID of subject
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
addpath([app_dir 'fieldtrip/']);
ft_defaults

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
if ischar(plot_scat); plot_scat = str2num(plot_scat); end

%% Prep variables
eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

SBJs = fn_load_SBJ_list(SBJ_id);

% Get condition info
[reg_lab, reg_names, reg_colors, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

% Load all ROI info
[roi_list, roi_colors, roi_field] = fn_roi_label_styles(roi_id);

% Set up electrode counts
reg_cnt  = zeros([numel(SBJs) numel(roi_list) numel(reg_lab)]);
elec_cnt = zeros([numel(SBJs) numel(roi_list)]);
out_cnt  = zeros(size(SBJs));
miss_cnt = zeros(size(SBJs));

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load ANOVA
    load([SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat']);
    
    %% Load ROI and GM/WM info
    load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    
%     % HACK!!!
%     cfgs = [];
%     if strcmp(SBJ,'CP24')
%         cfgs.channel = {'all','-RTO4'};
%     elseif strcmp(SBJ,'IR57')
%         cfgs.channel = {'all','-LAM4-5','-LAM5-6','-RIN8-9','-RIN9-10',...
%                         '-RSM1-2','-RSM2-3','-RSM3-4','-RTI9-10'};
%     elseif strcmp(SBJ,'IR68')
%         cfgs.channel = {'all','-LPC5-6','-LPC6-7','-LPC7-8'};
%     end
%     stat = ft_selectdata(cfgs,stat);

%     % Sort elecs by stat labels
%     if ~strcmp(SBJ,'IR66')  % HACK!!!
%         cfgs = []; cfgs.channel = stat.label;
%         elec = fn_select_elec(cfgs,elec);
%         elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
%     else
%         elec.roi = elec.(roi_id);
%     end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(beta.label)
        elec_ix = strcmp(beta.label{ch_ix},elec.label);
        if any(elec_ix)
            roi_ix = strcmp(elec.(roi_field){elec_ix},roi_list);
            if any(roi_ix)
                elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
                % Check for ANOVA group effects
                for reg_ix = 1:numel(reg_lab)
                    beta_reg_ix = strcmp(beta.feature,reg_lab{reg_ix});
                    if any(beta.qval(beta_reg_ix,ch_ix,:)<=st.alpha,3)
                        reg_cnt(sbj_ix,roi_ix,reg_ix) = reg_cnt(sbj_ix,roi_ix,reg_ix)+1;
                    end
                end
            else
                out_cnt(sbj_ix) = out_cnt(sbj_ix)+1;
                fprintf(2,'\t%s: %s in %s not in %s!\n',SBJ,beta.label{ch_ix},...
                    elec.(roi_field){elec_ix},roi_id);
            end
        else
            miss_cnt(sbj_ix) = miss_cnt(sbj_ix)+1;
            fprintf(2,'\t%s: %s in beta missing from elec!\n',SBJ,beta.label{ch_ix});
        end
    end
    clear SBJ SBJ_vars beta elec
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
if plot_scat
    scat_suffix = '_SBJscat';
else
    scat_suffix = '';
end
% Create and format the plot
fig_name = [SBJ_id '_HFA_errbar_' model_id '_' stat_id '_' atlas_id '_' roi_id scat_suffix];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% Compile Data in Plotting Format
scat_vars = cell(size(reg_lab));
bar_data  = zeros([numel(reg_lab) numel(roi_list)]);
sem_data  = zeros([numel(reg_lab) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
    for reg_ix = 1:numel(reg_lab)
        bar_data(reg_ix,roi_ix) = sum(reg_cnt(:,roi_ix,reg_ix))/sum(elec_cnt(:,roi_ix));
        sem_data(reg_ix,roi_ix) = nanstd(reg_cnt(:,roi_ix,reg_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
        scat_vars{reg_ix} = squeeze(reg_cnt(:,:,reg_ix));
        % Remove SBJs with no data in that ROI (make NaN)
        scat_vars{reg_ix}(elec_cnt==0) = nan;
    end
end

% Plot Activations by ROI
bars = cell(size(reg_lab));
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(roi_list));   %bar for activation, deactivation, condition
for reg_ix = 1:numel(reg_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    bars{reg_ix} = bar(bar_offsets+reg_ix,diag(bar_data(reg_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(bars{reg_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        
        % Plot standard error of the mean across subjects
        line([bar_offsets(roi_ix) bar_offsets(roi_ix)]+reg_ix,...
            [bar_data(reg_ix,roi_ix)+sem_data(reg_ix,roi_ix) bar_data(reg_ix,roi_ix)-sem_data(reg_ix,roi_ix)],...
            'Color','k','LineWidth',1.5);
        
        % Overlay scatter for individual SBJ
        if plot_scat
            has_elecs = squeeze(elec_cnt(:,roi_ix)~=0);
            scat = scatter(scat_offsets(has_elecs)+bar_offsets(roi_ix)+reg_ix,...
                scat_vars{reg_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),50,'k*');
        end
    end
end
leg_loc = 'northwest';
if plot_scat
    legend([bars{1},scat],roi_list{:},'Individuals','Location',leg_loc);
else
    legend(bars{1},roi_list{:},'Location',leg_loc);
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim       = [0.5 0.5+numel(reg_lab)];
ax.XTick      = 1:numel(reg_lab);
ax.XColor     = 'k';
ax.XTickLabel = reg_names;

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLabel.FontSize = 16;
ax.YLim            = [0 1];%ax.YLim(2)];%[0 0.6];%
ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Proportion of Electrodes Showing Significant Effects';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
set(gca,'FontSize',16);

%% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/errbar_ROI/'...
        model_id '/' stat_id '/' an_id '/' atlas_id '_' roi_id '/'];
    if ~exist(fig_dir,'dir'); [~,~] = mkdir(fig_dir); end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compute stats for differences
fprintf(2,'WARNING: t-test is likely not the right stat for this!\n');
pairs = nchoosek(1:numel(roi_list),2);
pvals = nan([numel(reg_lab) size(pairs,1)]);
for reg_ix = 1:numel(reg_lab)
    for p_ix = 1:size(pairs,1)
        [~, pvals(reg_ix,p_ix)] = ttest2(scat_vars{reg_ix}(:,pairs(p_ix,1)),...
                scat_vars{reg_ix}(:,pairs(p_ix,2)));
        if pvals(reg_ix,p_ix)<=0.05; sig_str = '*'; else sig_str = ''; end
        fprintf('%s%s: %s (%.2f) vs. %s (%.2f) p = %.5f\n',sig_str,...
            reg_lab{reg_ix},roi_list{pairs(p_ix,1)},bar_data(reg_ix,pairs(p_ix,1)),...
            roi_list{pairs(p_ix,2)},bar_data(reg_ix,pairs(p_ix,2)),pvals(reg_ix,p_ix));
    end
end

end
