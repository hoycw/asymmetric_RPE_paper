function BHV02d_RL_model_plot_grp_predictions_bar(SBJ_id,proc_id,model_id,plt_id,save_fig,varargin)
%% Plot model predictors across group averaged within condition
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   model_id [str] - ID of the model parameters to use
%   plt_id [str] - ID of the plotting parameters to use
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error/Apps/'];
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');     fig_vis = 'on'; end
if ~exist('fig_ftype','var');   fig_ftype = 'png'; end
if ischar(save_fig);            save_fig = str2num(save_fig); end

%% Load Data 
if ~strcmp(model_id,'ERPEs_DifFB'); error('run this for ERPEs_DifFB!'); end
proc_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
model_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m'];
eval(model_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
% [reg_lab, reg_names, reg_colors, reg_styles, reg_mrkrs] = fn_regressor_label_styles(mdl.model_lab);
if strcmp(mdl.model_lab,'RL3D')
    [reg_lab, reg_names, reg_colors, reg_styles] = fn_performanceRL_regressor_label_styles(mdl.model_lab);
else
    [reg_lab, reg_names, reg_colors, reg_styles, reg_mrkrs] = fn_regressor_label_styles(mdl.model_lab);
end
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(mdl.model_cond);
ez_idx = ~cellfun(@isempty,strfind(cond_lab,'Ez'));

%% Load Model and Comute Means
model = nan([numel(reg_lab) numel(cond_lab) numel(SBJs)]);
for s = 1:numel(SBJs)
    % Load RL Model
    load([root_dir 'PRJ_Error/data/' SBJs{s} '/03_events/' SBJs{s} '_bhv_' proc_id '_final.mat'],'bhv');
    tmp = load([root_dir 'PRJ_Error/data/' SBJs{s} '/06_models/' SBJs{s} '_model_' mdl.model_id '.mat']);
    
    % Average within condition
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        model(:,cond_ix,s) = nanmean(tmp.model(full_cond_idx==cond_ix,:),1);
    end
end

%% Compute Group Averages
% Variance is very small, so plot standard deviation instead of standard error of the mean
plot_means = nanmean(model,3);
plot_stds  = nan([numel(reg_lab) numel(cond_lab)]);
% plot_sems  = nan([numel(reg_lab) numel(cond_lab)]);
for reg_ix = 1:numel(reg_lab)
    plot_stds(reg_ix,:) = nanstd(model(reg_ix,:,:),[],3);
%     plot_sems(reg_ix,:) = nanstd(model(reg_ix,:,:),[],3)./sqrt(numel(SBJs))';
end

%% Plot Predictors
fig_name = [SBJ_id '_' mdl.model_id '_predictions_bar'];
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','OuterPosition',[0 0 0.5 0.5]);
ax = gca; hold on;

[~, sRPE_order] = sort(plot_means(strcmp(reg_lab,'sRPE'),:),'descend');

% Plot Activations by ROI
bars = cell(size(reg_lab));
bar_offsets  = linspace(-0.25,0.25,numel(cond_lab));   %bar for each condition
sbj_offsets = linspace(-0.02,0.02,numel(SBJs));
scat_sz = 25; scat_color = [0.3 0.3 0.3];
for reg_ix = 1:numel(reg_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    bars{reg_ix} = bar(bar_offsets+reg_ix,diag(plot_means(reg_ix,sRPE_order)),0.9,'stacked');
    for cond_ix = 1:numel(cond_lab)
        if contains(cond_lab(sRPE_order(cond_ix)),'Ez')
            set(bars{reg_ix}(cond_ix),'FaceColor',cond_colors{sRPE_order(cond_ix)},'EdgeColor','k');
        else
            set(bars{reg_ix}(cond_ix),'FaceColor',cond_colors{sRPE_order(cond_ix)},'EdgeColor','k','FaceAlpha',0.5);
        end
        
        % Plot SBJ mean datapoints
        scatter(sbj_offsets+bar_offsets(cond_ix)+reg_ix,model(reg_ix,sRPE_order(cond_ix),:),scat_sz,scat_color);
        
        % Plot standard error of the mean across subjects
        line([bar_offsets(cond_ix) bar_offsets(cond_ix)]+reg_ix,...
            [plot_means(reg_ix,sRPE_order(cond_ix))+plot_stds(reg_ix,sRPE_order(cond_ix)) ...
            plot_means(reg_ix,sRPE_order(cond_ix))-plot_stds(reg_ix,sRPE_order(cond_ix))],...
            'Color','k','LineWidth',1.5);
    end
end

% Plot parameters
set(gca,'XTick',1:numel(reg_lab));
set(gca,'XTickLabels',reg_names);
% xtickangle(plt.tick_angle);
xlim([0.5 numel(reg_lab)+0.5]);
ylim([-2 2]);

model_str = mdl.model_id;
title(model_str,'Interpreter','none');
legend(bars{1},cond_names{sRPE_order},'Location','southeast');
set(gca,'FontSize',16);

%% Save Figure
if save_fig
    % Create figure directory
    fig_dir = [root_dir 'PRJ_Error/results/model_predictions/' plt_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
