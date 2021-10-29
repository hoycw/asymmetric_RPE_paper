function SBJ07c_HFA_plot_stack_mean(SBJ, conditions, proc_id, an_id, actv_win,...
                                        plt_id, save_fig, varargin)
%% Plots single trial stack and condition averaged HFA per electrode with activation vs. baseline
%   sorts by condition, then by RT; scatter for events to show conditions
%   Optional: grand median instead of grand average
% INPUTS:
%   SBJ [str] - ID of subject
%   conditions [str] - - group of condition labels to segregate trials
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   actv_win [int] - length in ms of minimum difference from baseline to be active
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
%       plot_median [0/1] - binary flag to plot median ERP instead of mean
%           default: 0
% OUTPUTS:
%   saves figure

if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
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
        elseif strcmp(varargin{v},'plot_median')
            plot_median = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');   fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('atlas_id','var');  atlas_id = 'Dx'; end
if ~exist('plot_median','var'); plot_median = 0; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load Results
% Load RTs
load([SBJ_vars.dirs.events SBJ '_bhv_' proc_id '_final.mat'],'bhv');

% Load data
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'.mat');
load(hfa_fname,'hfa');
actv_fname = strcat(hfa_fname(1:end-4),'_actv_mn',num2str(actv_win),'.mat');
load(actv_fname);

% Load ROI and GM/WM info
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat'];
try
    load(elec_fname);
catch
    fprintf('Final elec file not available, no ROI info in title: %s\n',elec_fname);
end

%% Select Data
% Select conditions
[cond_lab, cond_names, cond_colors, cond_styles, cond_mrkrs] = fn_condition_label_styles(conditions);
full_cond_idx = fn_condition_index(cond_lab,bhv);
bhv = fn_select_bhv(bhv, full_cond_idx);
cond_idx = fn_condition_index(cond_lab, bhv);

% Sort trials by condition, RT, then trial number
cond_mat   = horzcat(cond_idx,bhv.rt,[1:numel(bhv.trl_n)]');
cond_mat   = sortrows(cond_mat,[1 2]);
cond_edges = find(diff(cond_mat(:,1)));

% Select channels in elec file
cfg_trim = [];
if exist('elec_lab','var')
    cfg_trim.channel = elec_lab;
    if exist('elec','var')
        elec = fn_select_elec(cfg_trim,elec);
    end
end

% Select channels, trials, and epochs in data
cfg_trim.trials  = find(full_cond_idx);
cfg_trim.latency = plt.plt_lim;
hfa = ft_selectdata(cfg_trim,hfa);

% Get event times for plotting
[evnt_times] = fn_get_evnt_times(an.evnt_lab,plt.evnt_lab,bhv);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error/results/HFA/' SBJ '/stack_mn_' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir'); [~] = mkdir(fig_dir); end
sig_ln_dir = [fig_dir 'sig_ch/'];
if ~exist(sig_ln_dir,'dir'); [~] = mkdir(sig_ln_dir); end

% Create a figure for each channel
for ch_ix = 1:numel(hfa.label)
    sig_flag = 0;
    %% Compute Plotting Data
    % Average HFA per condition
    means = NaN([numel(cond_lab) numel(hfa.time)]);
    sems  = NaN([numel(cond_lab) numel(hfa.time)]);
    for cond_ix = 1:numel(cond_lab)
        if plot_median
            means(cond_ix,:) = squeeze(median(hfa.powspctrm(cond_idx==cond_ix,ch_ix,1,:),1));
        else
            means(cond_ix,:) = squeeze(mean(hfa.powspctrm(cond_idx==cond_ix,ch_ix,1,:),1));
        end
        sems(cond_ix,:) = squeeze(std(hfa.powspctrm(cond_idx==cond_ix,ch_ix,1,:),[],1))./sqrt(sum(cond_idx==cond_ix))';
    end
    
    %% Create Plot and Parameters
    fig_name = [SBJ '_' conditions '_stack_' hfa.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 1],'Visible',fig_vis);
    
    % Get color limits
    clims = NaN([1 2]);
    clims(1) = prctile(reshape(hfa.powspctrm(:,ch_ix,:,:),...
        [1 size(cond_mat,1)*numel(hfa.time)]),plt.clim_perc(1));
    clims(2) = prctile(reshape(hfa.powspctrm(:,ch_ix,:,:),...
        [1 size(cond_mat,1)*numel(hfa.time)]),plt.clim_perc(2));
    clims = [min(clims(1)) max(clims(2))];
    
    %% Plot Single Trial Stack
    subplot(3,1,1:2); hold on;
    ax = gca;
    
    % Plot Single Trials Per Condition
    imagesc(hfa.time,1:size(cond_mat,1),squeeze(hfa.powspctrm(cond_mat(:,3),ch_ix,:,:)));
    set(gca,'YDir','normal');
    
    % Plot events: stim, target, feedback onset, feedback offset
    evnt_lines = gobjects(size(evnt_times));
    for evnt_ix = 1:numel(evnt_times)
        evnt_lines(evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt.evnt_width,'Color','k','LineStyle',plt.evnt_styles{evnt_ix});
    end
    
    % Scatter plot of RTs or feedback onset to mark conditions
    trl_thin_factor = 15;%floor(size(cond_mat,1)/100);    % for feedback plotting, thin out the RTs to see them
    scat = gobjects(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        idx = cond_mat(:,1)==cond_ix;
        if any(strcmp(plt.evnt_lab,'R'))
            scat(cond_ix) = scatter(cond_mat(idx,2),find(idx),plt.evnt_mrkr_sz,...
                'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k',...
                'Marker',cond_mrkrs{cond_ix});
        elseif any(strcmp(plt.evnt_lab,'F'))
            trl_ix = find(idx);
            trl_ix = trl_ix(1:trl_thin_factor:end);
            scat(cond_ix) = scatter(zeros(size(trl_ix)),trl_ix,plt.evnt_mrkr_sz,...
                'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k',...
                'Marker',cond_mrkrs{cond_ix});
        end
    end
    ylim([1 size(cond_mat,1)]);
    
    % Plot condition dividers
    for div_ix = 1:numel(cond_edges)
        line(plt.plt_lim, [cond_edges(div_ix) cond_edges(div_ix)],...
            'LineWidth',plt.evnt_width,'Color','k','LineStyle','--');
    end
    
    % Plotting parameters
    title_str = hfa.label{ch_ix};
    if exist('elec','var')
        title_str  = [title_str ' (' elec.ROI{strcmp(hfa.label{ch_ix},elec.label)} ')'];
    end
    ax.Title.String  = title_str;
    ax.XLim          = plt.plt_lim;
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
%     ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Trials';
    cbar = colorbar;
    caxis(clims);
    if plt.legend
        legend(scat,cond_names{:}, 'Location',plt.legend_loc);
%         if ~isempty(scat)
%         else
%             legend(evnt_lines, plt.evnt_lab, 'Location',plt.legend_loc);
%         end
%         legend(event_lines,event_lab,'Location',plt.legend_loc);
    end
    set(gca,'FontSize',16);
    
    %% Plot Condition Averaged HFA
    subplot(3,1,3); hold on;
    ax2 = gca;
    
    % Plot HFA Means (and variance)
    cond_lines = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(plt.evnt_lab) 1]);
    main_line_ix = 0;
    for cond_ix = 1:numel(cond_lab)
        main_line_ix = main_line_ix + 1;
        cond_lines{cond_ix} = shadedErrorBar(hfa.time, means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(main_line_ix) = cond_lines{cond_ix}.mainLine;
    end
    
    % Plot Significance for Activation vs. Baseline
    ylims = ylim;%[-0.3 1.3];
    if actv.actv_ch(ch_ix)
        sig_flag = 1;
        for sig_ix = 1:size(actv.actv_epochs{ch_ix},1)
            % Plot rectangular patch over epoch
            patch([actv.actv_epochs{ch_ix}(sig_ix,1) actv.actv_epochs{ch_ix}(sig_ix,1) ...
                actv.actv_epochs{ch_ix}(sig_ix,2) actv.actv_epochs{ch_ix}(sig_ix,2)],...
                [ylims(1) ylims(2) ylims(2) ylims(1)],...
                plt.sig_color,'FaceAlpha',plt.sig_alpha);
        end
    end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_line_ix = main_line_ix + 1;
        main_lines(main_line_ix) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],ylims,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end    
    
    % Axes and Labels
    ax2.YLabel.String = 'HFA Power (z)';
    ax2.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax2.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax2.XLabel.String = 'Time (s)';
    ax2.Title.String  = 'Condition Averaged HFA';
    if plt.legend
        legend(main_lines,[cond_lab plt.evnt_lab],'Location',plt.legend_loc);
    end
    set(gca,'FontSize',16);
    ax2.YLim = ylims;
    
    %% Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        % Ensure vector graphics if saving
        if any(strcmp(fig_ftype,{'svg','eps'})); set(gcf, 'Renderer', 'painters'); end
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
