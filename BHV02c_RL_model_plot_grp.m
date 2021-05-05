function BHV02c_RL_model_plot_grp(SBJ_id,proc_id,model_id, varargin)
%% Plot all SBJ RL model fits overlapping
%   Rereferenced to the midpoint tolerance
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
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

%% Load Parameters
proc_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
model_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m'];
eval(model_vars_cmd);

% Load SBJ list
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(mdl.model_cond);

% Initialize Plotting Variables
sig_step = 0.001;               % tolerance step size for plotting model fit
sig_x    = [0:sig_step:0.4];
sig_y    = nan([numel(sig_x) numel(SBJs)]);

%% Load Data
for s = 1: numel(SBJs)
    % Load Subject Specific Data
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
    load([SBJ_vars.dirs.events SBJ '_bhv_' proc_id '_final.mat'],'bhv');

    % Select Conditions of Interest
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    orig_n_trials = numel(bhv.trl_n);
    bhv = fn_select_bhv(bhv, full_cond_idx);
    fprintf('\t%s: Loaded %d trials, kept %d for modeling...\n',SBJ,orig_n_trials,numel(bhv.trl_n));
    
    % Recalculate Model Fit
    sig_y(:,s) = betas(1) + (sig_x * betas(2));
    sig_y(:,s) = 1 ./ (1+exp(-sig_y(:,s)));
    
    clear SBJ_vars
end

%% Plot All SBJ Level Sigmoids
fig_name = ['GRP_BHV_acc_' mdl.model_id '_pWin'];
fig = figure('Name',fig_name,'Visible',fig_vis);
hold on;

% Plot model fits
for s = 1:numel(SBJs)
    line(sig_x,sig_y(:,s),'Color', 'k');
end

% Figure Parameters
title('GRP Level Accuracy vs. Tolerance');
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'YLim',[0 1]);
set(gca,'FontSize',14);

%% Save Figures
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
    saveas(fig,fig_fname);
end

end
