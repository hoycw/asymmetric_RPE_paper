function BHV02b_RL_model_plot(SBJ,proc_id,model_id,varargin)
%% Plot single SBJ behavior with RL model fit
%   Scatter of single trial tolerance and outcomes
%   Scatter of block accuracy and average tolerance
%   Sigmoid from logistic regression fit on top
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   model_id [str] - ID of the model parameters to use
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       save_fig [0/1] - binary flag to save figure; default = 1
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error/Apps/'];
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig') && ischar(varargin{v+1})
            save_fig = varargin{v+1};
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
if ~exist('save_fig','var');    save_fig = 1; end

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
model_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m'];
eval(model_vars_cmd);

% Plotting parameters
trl_sz = 25;        % scatter size for single trial
mean_sz = 100;      % scatter size for blocks
sig_step = 0.001;   % tolerance step size for plotting model fit

% Get model and condition parameters
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(mdl.model_cond);

%% Load and Select Behavior
% Load data (should include betas from logistic regression)
load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);

load([SBJ_vars.dirs.events SBJ '_bhv_' proc_id '_final.mat'],'bhv');

% Select Trails based on Conditions of Interest
full_cond_idx = fn_condition_index(cond_lab, bhv);
orig_n_trials = numel(bhv.trl_n);
bhv = fn_select_bhv(bhv, full_cond_idx);
n_trials = numel(bhv.trl_n);
fprintf('\t%s: Loaded %d trials, kept %d for modeling...\n',SBJ,orig_n_trials,n_trials);

%% Compute Mean Accuracy and Tolerance per Block
% % Adjust block numbers if cut off in middle of task
% if strcmp(SBJ,'EEG12')
%     blk5_starts = find(bhv.blk==5 & bhv.blk_trl_n==1);
%     for trl_ix = 1:blk5_starts(2)-1
%         bhv.blk(trl_ix) = bhv.blk(trl_ix)-4;
%     end
% end

% Compute block level accuracy and tolerance
run_ids = unique(bhv.run);
blk_ids = unique(bhv.blk);
blk_ez  = false([numel(run_ids) numel(blk_ids)]);
blk_tol = nan([numel(run_ids) numel(blk_ids)]);
blk_acc = nan([numel(run_ids) numel(blk_ids)]);
for r_ix = 1:numel(run_ids)
    for b_ix = 1:numel(blk_ids)
        blk_tol(r_ix,b_ix) = mean(bhv.tol(bhv.run==run_ids(r_ix) & bhv.blk==blk_ids(b_ix)));
        blk_acc(r_ix,b_ix) = mean(bhv.hit(bhv.run==run_ids(r_ix) & bhv.blk==blk_ids(b_ix)));
        if strcmp(unique(bhv.cond(bhv.run==run_ids(r_ix) & bhv.blk==blk_ids(b_ix))),'easy')
            blk_ez(r_ix,b_ix) = true;
        end
    end
end

%% Plot Tolerance vs. Outcome with Model Overlay
fig_name = [SBJ '_BHV_acc_' mdl.model_id '_pWin'];
figure('Name',fig_name,'Visible',fig_vis);
hold on;

% Plot Trial and Block Behavior: Easy
ez_trl_idx = strcmp(bhv.cond,'easy');
ez_trl = scatter(bhv.tol(ez_trl_idx), bhv.hit(ez_trl_idx), trl_sz,'k','filled');
ez_blk = scatter(blk_tol(blk_ez),blk_acc(blk_ez), mean_sz,'k','filled');

% Plot Trial and Block Behavior: Hard
hd_trl = scatter(bhv.tol(~ez_trl_idx), bhv.hit(~ez_trl_idx), trl_sz,'k','d');
hd_blk = scatter(blk_tol(~blk_ez),blk_acc(~blk_ez), mean_sz,'k','d');

% Reconstruct and plot Model Fit
sig_x = [0:sig_step:0.4];
sig_y = betas(1) + (sig_x * betas(2));
sig_y = 1 ./ (1+exp(-sig_y));
fit_line = line(sig_x,sig_y,'Color','k');

% Figure Parameters
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'YLim',[0 1]);
title(SBJ);
ez_leg = ['Easy (n=' num2str(sum(ez_trl_idx)) '; mean = ' ...
    num2str(mean(bhv.hit(ez_trl_idx))*100,'%.1f') '%)'];
hd_leg = ['Hard (n=' num2str(sum(~ez_trl_idx)) '; mean = ' ...
    num2str(mean(bhv.hit(~ez_trl_idx))*100,'%.1f') '%)'];
legend([ez_trl, hd_trl, fit_line],{ez_leg,hd_leg,'Model Fit'},'Location','southeast');
set(gca,'FontSize',14);

%% Save Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error/results/BHV/model_fits/' mdl.model_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(gcf,fig_fname);
end

end
