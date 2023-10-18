function SBJ08c_HFA_plot_grp_mGLM_ROI_hist(SBJ_id,proc_id,an_id,model_id,stat_id,...
                                    roi_id,save_fig,varargin)
% Plots histograms of mGLM max beta by ROI
%       COMMENTS NOT YET ADJUSTED
%   gROI version: one plot; subplots per gROI (lines still colored by ROI)
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   model_id [str] - ID of the model to load
%   stat_id [str] - ID of statistical analysis
%   an_id [str] - analysis ID for S-locked preprocessing, filtering, etc.
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
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
% eval(['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);

SBJs = fn_load_SBJ_list(SBJ_id);

% Get condition info
[reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

% Load all ROI info
[roi_list, roi_colors, roi_field] = fn_roi_label_styles(roi_id);
if ~strcmp(roi_field,'gROI'); error('not ready for non gROI'); end

%% Load Data
max_beta  = cell(size(SBJs));
beta_sig  = cell(size(SBJs));
beta_sig_tp  = cell(size(SBJs));    % specifically that time point significant?
beta_most_sig  = cell(size(SBJs));    % is the max beta time point the most sig time point?
roi_ch    = cell(size(SBJs));
onsets    = cell([numel(reg_lab) numel(SBJs)]);
roi_cnt   = zeros(size(roi_list));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([SBJ_vars.dirs.recon SBJs{s} '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    load([SBJ_vars.dirs.stats SBJs{s} '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat']);
    
    max_beta{s} = nan([numel(beta.label) numel(reg_lab)]);
    beta_sig{s} = false([numel(beta.label) numel(reg_lab)]);
    beta_sig_tp{s} = false([numel(beta.label) numel(reg_lab)]);
    beta_most_sig{s} = false([numel(beta.label) numel(reg_lab)]);
    roi_ch{s}   = zeros(size(beta.label));
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
                    [~,max_ix] = max(abs(beta.trial(beta_reg_ix,ch_ix,:)));
                    max_beta{s}(ch_ix,reg_ix) = beta.trial(beta_reg_ix,ch_ix,max_ix);
                    if any(beta.qval(beta_reg_ix,ch_ix,:)<=0.05,3)
                        beta_sig{s}(ch_ix,reg_ix) = true;
%                         onsets{reg_ix,roi_ix} = [onsets{reg_ix,roi_ix}...
%                             min(time_vec(beta.qval(beta_reg_ix,ch_ix,:)<=st.alpha))];
                    end
                    if beta.qval(beta_reg_ix,ch_ix,max_ix)<=0.05
                        beta_sig_tp{s}(ch_ix,reg_ix) = true;
                    end
                    [~,min_ix] = min(beta.pval(beta_reg_ix,ch_ix,:));
                    if max_ix==min_ix
                        beta_most_sig{s}(ch_ix,reg_ix) = true;
                    end
                end
            end
        else
            fprintf(2,'\t%s: %s in beta missing from elec!\n',SBJ,beta.label{ch_ix});
        end
    end
    
    clear beta elec SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine across SBJs
grp_max_beta      = cat(1,max_beta{:});
grp_beta_sig      = cat(1,beta_sig{:});
grp_beta_sig_tp   = cat(1,beta_sig_tp{:});
grp_beta_most_sig = cat(1,beta_most_sig{:});
grp_roi_ch        = cat(1,roi_ch{:});
nan_idx = isnan(grp_max_beta(:,1));

% Check correspondence between max beta and significance
for reg_ix = 1:numel(reg_lab)
    fprintf('%s: %d / %d max betas are significant\n',reg_lab{reg_ix},sum(grp_beta_sig(~nan_idx,reg_ix)),sum(~nan_idx));
    fprintf('%s: %d / %d max betas are significant at that time point\n',reg_lab{reg_ix},sum(grp_beta_sig_tp(~nan_idx,reg_ix)),sum(~nan_idx));
    fprintf('%s: %d / %d max betas occur at most significant time point\n',reg_lab{reg_ix},sum(grp_beta_most_sig(~nan_idx,reg_ix)),sum(~nan_idx));
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/GLM_beta_hist/' model_id '/' stat_id '/' an_id '/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

for roi_ix = 1:numel(roi_list)
    % Create a figure for each condition (all gROIs in one subplot)
    fig_name = [SBJ_id '_GLM_beta_hist_' model_id '_' stat_id '_' atlas_id '_' roi_list{roi_ix}];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis); %twice as wide for the double plot
    
    % Set xlims per ROI
    n_bins = 100;
    all_betas = grp_max_beta(~nan_idx,:);
    lim_edge = max(abs(all_betas(:)));
    hist_bins = linspace(-lim_edge,lim_edge,n_bins);
    
    for reg_ix = 1:numel(reg_lab)
        ax = subplot(numel(reg_lab),1,reg_ix);%subplot_ix);
        hold on;
        % Plot all max betas
        histogram(grp_max_beta((grp_roi_ch==roi_ix & ~nan_idx),reg_ix),hist_bins,'FaceColor','k','FaceAlpha',0.1);
        % Plot significant max betas
        histogram(grp_max_beta(grp_roi_ch==roi_ix & grp_beta_sig(:,reg_ix)==1,reg_ix),hist_bins,'FaceColor',roi_colors{roi_ix},'FaceAlpha',0.9);
        
        
        % Plotting parameters
        ax.Title.String  = [reg_names{reg_ix} ' (n = ' num2str(sum(grp_beta_sig(grp_roi_ch==roi_ix,reg_ix))) '/' num2str(roi_cnt(roi_ix)) ')'];
        ax.XLabel.String = 'Max Beta Coefficient';
        set(gca,'FontSize',16);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Plot scatter of shared significant coefficients
% Create a figure for each condition (all gROIs in one subplot)
fig_name = [SBJ_id '_GLM_beta_scatter_' model_id '_' stat_id '_' atlas_id '_' roi_id];
f = figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.4 1],'Visible',fig_vis); %twice as wide for the double plot

pairs = nchoosek(1:numel(reg_lab),2);

roi_mrkrs = {'*','o','d','x'};
for p_ix = 1:size(pairs,1)
    ax = subplot(size(pairs,1),1,p_ix);
    hold on;
    dbl_sig_idx = [grp_beta_sig(:,pairs(p_ix,1))==1 & grp_beta_sig(:,pairs(p_ix,2))==1];
    for roi_ix = 1:numel(roi_list)
        scatter(grp_max_beta(dbl_sig_idx & grp_roi_ch==roi_ix,pairs(p_ix,1)),...
                grp_max_beta(dbl_sig_idx & grp_roi_ch==roi_ix,pairs(p_ix,2)),...
                'Marker',roi_mrkrs{roi_ix},'MarkerEdgeColor',roi_colors{roi_ix},...
                'LineWidth',2,'SizeData',100);
    end
    line([-lim_edge lim_edge],[-lim_edge lim_edge],'Color','k','LineStyle','--');
    line([-lim_edge lim_edge],[lim_edge -lim_edge],'Color','k','LineStyle','--');
    line([0 0],[-lim_edge lim_edge],'Color','k','LineStyle','-');
    line([-lim_edge lim_edge],[0 0],'Color','k','LineStyle','-');
    ax.XLim = [-lim_edge lim_edge];
    ax.YLim = [-lim_edge lim_edge];
    ax.XLabel.String = reg_names{pairs(p_ix,1)};
    ax.YLabel.String = reg_names{pairs(p_ix,2)};
    ax.FontSize = 16;
end

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
