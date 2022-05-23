function SBJ08h_HFA_plot_crRT_LME(SBJ, proc_id, an_id, model_id, stat_id, plt_id, save_fig, varargin)
%% Plots LME beta time series per electrode
% INPUTS:
%   SBJ [str] - ID of subject
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   model_id [str] - ID of the model used in GLM
%   stat_id [str] - ID of the statistical parameters to extract from HFA
%   plt_id [str] - ID of the plotting parameters to use
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       atlas_id [str] - ID of the atlas for elec ROI assignments
%           default: 'Dx'
%       elec_lab [cell array] - list of electrode labels to plot only subset
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
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
        elseif strcmp(varargin{v},'elec_lab') && iscell(varargin{v+1})
            elec_lab = varargin{v+1};
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

%% Load Results
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);

[reg_lab, reg_names, reg_colors, reg_styles, ~] = fn_regressor_label_styles(mdl.model_lab);

% Load RTs
load([SBJ_vars.dirs.events SBJ '_bhv_' proc_id '_final.mat'],'bhv');
stats_fname = [root_dir 'PRJ_Error/data/GRP/stats/' model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname);

% Load ROI and GM/WM info
elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_' atlas_id '_final.mat'];
load(elec_fname);

% Get event times for plotting
[evnt_times] = fn_get_evnt_times(an.evnt_lab,plt.evnt_lab,bhv);

%% Select SBJ (or electrode) specifci data
if exist('elec_lab','var')
    beta.coefs   = zeros(numel(elec_lab),size(beta_chan.coefs{1},2),size(beta_chan.coefs{1},3));
    beta.qvals   = nan(numel(elec_lab),size(beta_chan.coefs{1},2),size(beta_chan.coefs{1},3));
    beta.label   = elec_lab;
    beta.time    = beta_chan.time;
    beta.feature = beta_chan.feature;
    beta.chancat = cell(size(elec_lab));
    for e_ix = 1:numel(elec_lab)
        % Get data
        for roi_ix = 1:numel(beta_chan.label)
            if any(contains(beta_chan.chan_label{roi_ix},[SBJ ' ' elec_lab{e_ix}]))
                ch_ix = find(strcmp(beta_chan.chan_label{roi_ix},[SBJ ' ' elec_lab{e_ix}]));
                beta.coefs(e_ix,:,:) = beta_chan.coefs{roi_ix}(ch_ix,:,:);
                beta.qvals(e_ix,:,:) = beta_chan.qvals{roi_ix}(ch_ix,:,:);
                % Get category
                for cat_ix = 1:numel(beta_chan.chancat_label)
                    cat_elecs = beta_chan.chancat_ix{roi_ix}{cat_ix};
                    if any(cat_elecs==ch_ix)
                        beta.chancat{e_ix}   = beta_chan.chancat_label{cat_ix};
                    end
                end
            end
        end
    end
else
    error('have not set up whole SBJ LME beta time series plotting yet');
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/' SBJ '/' model_id '/' stat_id '/' an_id '/' plt_id '/'];
if ~exist(fig_dir,'dir'); [~] = mkdir(fig_dir); end
sig_ln_dir = [fig_dir 'sig_ch/'];
if ~exist(sig_ln_dir,'dir'); [~] = mkdir(sig_ln_dir); end

% Find plot limits
max_beta = max(beta.coefs(:));
min_beta = min(beta.coefs(:));
ylim_fudge = (max_beta-min_beta)*plt.ylim_fudge;
ylims  = [min_beta-ylim_fudge max_beta+ylim_fudge];

% Create a figure for each channel
for ch_ix = 1:numel(beta.label)
    if ~exist('elec_lab','var') || (exist('elec_lab','var') && any(strcmp(beta.label{ch_ix},elec_lab)))
        sig_flag = 0;
        fig_name = [SBJ '_' model_id '_' stat_id '_' beta.label{ch_ix}];
        figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 1 0.8],'Visible',fig_vis);
        
        ax = gca; hold on;
        main_lines = gobjects(size(reg_lab));
        % Plot var_exp
        for reg_ix = 1:numel(reg_lab)
            beta_reg_ix = find(strcmp(beta.feature,reg_lab{reg_ix}))+1;
            main_lines(reg_ix) = plot(beta.time,squeeze(beta.coefs(ch_ix,beta_reg_ix,:)),...
                'Color',reg_colors{reg_ix},'LineStyle',reg_styles{reg_ix});
        end
        
        % Plot events: stim, target, feedback onset, feedback offset
        evnt_lines = gobjects(size(plt.evnt_lab));
        for evnt_ix = 1:numel(evnt_times)
            evnt_lines(evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylims,...
                'LineWidth',plt.evnt_width,'Color','k','LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Plot significant time periods
        for reg_ix = 1:numel(reg_lab)
            beta_reg_ix = find(strcmp(beta.feature,reg_lab{reg_ix}))+1;
            if any(squeeze(beta.qvals(ch_ix,beta_reg_ix,:))<=st.alpha)
                sig_flag = 1;
                % Find significant periods
                sig_chunks = fn_find_chunks(squeeze(beta.qvals(ch_ix,beta_reg_ix,:))<=st.alpha);
                sig_chunks(squeeze(beta.qvals(ch_ix,beta_reg_ix,sig_chunks(:,1)))>st.alpha,:) = [];
                fprintf('%s %s -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                    beta.label{ch_ix},reg_lab{reg_ix},size(sig_chunks,1));
                
                % Plot Significance
                for sig_ix = 1:size(sig_chunks,1)
                    if strcmp(plt.sig_type,'bold')
                        if diff(sig_chunks(sig_ix,:))==0    % single windows
                            scatter(beta.time(sig_chunks(sig_ix,1)),squeeze(beta.coefs(ch_ix,beta_reg_ix,sig_chunks(sig_ix,1))),...
                                plt.sig_scat_size,reg_colors{reg_ix},plt.sig_scat_mrkr,'filled');
                        else
                            line(beta.time(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                                squeeze(beta.coefs(ch_ix,beta_reg_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                                'Color',reg_colors{reg_ix},'LineStyle',plt.sig_style,...
                                'LineWidth',plt.sig_width);
                        end
                    elseif strcmp(plt.sig_type,'patch')
                        error('sig_type = patch needs sig_times variable!');
                        sig_times = win_lim(sig_chunks(sig_ix,:),:);
                        patch([sig_times(1,1) sig_times(1,1) sig_times(2,2) sig_times(2,2)], ...
                            [ylims(1) ylims(2) ylims(2) ylims(1)],...
                            reg_colors{reg_ix},'FaceAlpha',plt.sig_alpha);
                    end
                end
            else
                %             fprintf('%s %s -- NO SIGNIFICANT CLUSTERS FOUND...\n',...
                %                 beta.label{ch_ix},reg_lab{reg_ix});
            end
        end
        
        % Plotting parameters
        elec_ix = find(strcmp(beta.label{ch_ix},elec.label));
        if ~isempty(elec_ix)
            ax.Title.String  = strcat(beta.label{ch_ix}, ' (', elec.ROI{elec_ix},...
                '-', beta.chancat{ch_ix} ,')');
        else
            ax.Title.String  = strcat(beta.label{ch_ix});
        end
        ax.Box           = 'off';
        ax.YLim          = ylims;
        ax.YLabel.String = 'Model Coefficient';
        ax.XLim          = plt.plt_lim;
        ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
        ax.XLabel.String = 'Time (s)';
        legend([main_lines evnt_lines],[reg_names plt.evnt_lab],'Location',plt.legend_loc);
        set(ax,'FontSize',16);
        
        % Save figure
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
            
            % Symbolic link for significant plots
            if sig_flag
                cd(sig_ln_dir);
                link_cmd = ['ln -s ../' fig_name '.' fig_ftype ' .'];
                system(link_cmd);
            end
        end
    end
end

end
