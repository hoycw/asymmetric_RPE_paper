function SBJ06g_ERP_plot_grp_GLM_onsets_ROI(SBJ_id,proc_id,an_id,model_id,stat_id,roi_id,...
                                            plt_id,save_fig,varargin)%gm_thresh,z_thresh,
% Load ERP mGLM results to plot onsets of signifcant effects per ROI
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   plt.grp_metric [str] - {'avg','mdn','all'}
%       mean/median will compute that metric within each SBJ (variance is across SBJs)
%       all- all electrode onsets are aggregated as if from the same SBJ

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

%% Prep variables
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);

SBJs = fn_load_SBJ_list(SBJ_id);

% Get condition info
[reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

% Load all ROI info
[roi_list, roi_colors, roi_field] = fn_roi_label_styles(roi_id);
if strcmp(plt.grp_metric,'all')
    SBJ_colors = distinguishable_colors(numel(SBJs));
end

%% Load Results
% Set up onset counts
onsets  = cell([numel(reg_lab) numel(roi_list) numel(SBJs)]);
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Load ROI and GM/WM info
    load([SBJ_vars.dirs.recon SBJs{s} '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
        
%     % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
%     gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
        
    %% Aggregate results per ROI
    % Load data
    load([SBJ_vars.dirs.stats SBJs{s} '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat'],'beta');
    if s==1; time_vec = beta.time; end
    for ch_ix = 1:numel(beta.label)
        elec_ix = strcmp(beta.label{ch_ix},elec.label);
        if any(elec_ix) %&& beta.max_hfa_z(ch_ix)>=z_thresh% && gm_bin(ch_ix)
            roi_ix = find(strcmp(roi_list,elec.(roi_field){elec_ix}));
            if ~isempty(roi_ix)
                % Consolidate to binary sig/non-sig
                for reg_ix = 1:numel(reg_lab)
                    beta_reg_ix = strcmp(beta.feature,reg_lab{reg_ix});
                    if any(beta.qval(beta_reg_ix,ch_ix,:)<=0.05,3)
                        onsets{reg_ix,roi_ix,s} = [onsets{reg_ix,roi_ix,s}...
                            min(time_vec(beta.qval(beta_reg_ix,ch_ix,:)<=st.alpha))];
                    end
                end
            end
        end
    end
    clear SBJ SBJ_vars elec beta
end

%% Aggregate/Process onsets per gROI
% Format as struct to fit violinplot
plot_onsets    = cell([numel(reg_lab) 1]);
plot_onset_sbj = cell([numel(reg_lab) 1]);
good_roi_ix    = cell([numel(reg_lab) 1]);
for reg_ix = 1:numel(reg_lab)
    for roi_ix = 1:numel(roi_list)
        if ~isempty([onsets{reg_ix,roi_ix,:}])
            if strcmp(plt.grp_metric,'all')
                plot_onsets{reg_ix}.(roi_list{roi_ix}) = [onsets{reg_ix,roi_ix,:}]';
                plot_onset_sbj{reg_ix}.(roi_list{roi_ix}) = [];
            end
            for s = 1:numel(SBJs)
                % Aggregate onsets per ROI within each SBJ
                if strcmp(plt.grp_metric,'all')
                    plot_onset_sbj{reg_ix}.(roi_list{roi_ix}) = ...
                        [plot_onset_sbj{reg_ix}.(roi_list{roi_ix}); repmat(s,size(onsets{reg_ix,roi_ix,s}))'];
                elseif strcmp(plt.grp_metric,'mdn')
                    plot_onsets{reg_ix}.(roi_list{roi_ix}) = ...
                        [plot_onsets{reg_ix}.(roi_list{roi_ix}); nanmedian(onsets{reg_ix,roi_ix,s})];
                elseif strcmp(plt.grp_metric,'avg')
                    plot_onsets{reg_ix}.(roi_list{roi_ix}) = ...
                        [plot_onsets{reg_ix}.(roi_list{roi_ix}); nanmean(onsets{reg_ix,roi_ix,s})];
                else
                    error(['Unknown plt.grp_metric: ' plt.grp_metric]);
                end
            end
        else
            plot_onsets{reg_ix} = struct();
            fprintf(2,'\tNo significant elecs for %s in %s\n',reg_names{reg_ix},roi_list{roi_ix});
        end
    end
    good_roi_ix{reg_ix} = find(contains(roi_list,fieldnames(plot_onsets{reg_ix})));
end

%% Pair-Wise Statistics
pairs = nchoosek(1:numel(roi_list),2);
pvals = nan([numel(reg_lab) size(pairs,1)]);
for reg_ix = 1:numel(reg_lab)
    for p_ix = 1:size(pairs,1)
        if isfield(plot_onsets{reg_ix},roi_list{pairs(p_ix,1)}) && isfield(plot_onsets{reg_ix},roi_list{pairs(p_ix,1)})
            [~, pvals(reg_ix,p_ix)] = ttest2(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,1)}),...
                plot_onsets{reg_ix}.(roi_list{pairs(p_ix,2)}));
            if pvals(reg_ix,p_ix)<=0.05; sig_str = '*'; else sig_str = ''; end
            fprintf('\t%s%s: %s (nE=%d; nSBJ=%d; %.3f +/- %.4f) vs. %s (nE=%d; nSBJ=%d; %.3f +/- %.4f) p = %.5f\n',...
                sig_str,reg_names{reg_ix},roi_list{pairs(p_ix,1)},...
                numel(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,1)})),...
                numel(unique(plot_onset_sbj{reg_ix}.(roi_list{pairs(p_ix,1)}))),...
                mean(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,1)})),...
                std(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,1)})), roi_list{pairs(p_ix,2)}, ...
                numel(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,2)})),...
                numel(unique(plot_onset_sbj{reg_ix}.(roi_list{pairs(p_ix,2)}))),...
                mean(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,2)})), std(plot_onsets{reg_ix}.(roi_list{pairs(p_ix,2)})),...
                pvals(reg_ix,p_ix));
        else
            fprintf('\t%s: pair %s vs. %s not possible\n',reg_names{reg_ix},...
                roi_list{pairs(p_ix,1)},roi_list{pairs(p_ix,2)});
        end
    end
end

%% Plot GROI Results
for reg_ix = 1:numel(reg_lab)
    % Create and format the plot
    fig_name = [SBJ_id '_ERP_onsets_' model_id '_' stat_id '_' roi_id '_' reg_lab{reg_ix}];
%         '_GM' num2str(gm_thresh) '_z' num2str(z_thresh) '_normRTout'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
    
    violins = violinplot(plot_onsets{reg_ix});%,'ViolinAlpha',0.3);
                        
    % Adjust plot propeties
    for roi_ix = good_roi_ix{reg_ix}
        violin_ix = find(good_roi_ix{reg_ix}==roi_ix);
        if isfield(plt,'violin_scat_colors') && strcmp(plt.violin_scat_colors,'SBJ')
            violins(violin_ix).ViolinColor = [0.8 0.8 0.8];
            violins(violin_ix).BoxPlot.FaceColor = roi_colors{roi_ix};
            violins(violin_ix).EdgeColor = roi_colors{roi_ix};
            
            % Change scatter colors wihtin violin to mark SBJ
            scat_colors = zeros([numel(violins(violin_ix).ScatterPlot.XData) 3]);
            for s = 1:numel(SBJs)
                scat_colors(plot_onset_sbj{reg_ix}.(roi_list{roi_ix})==s,:) = repmat(SBJ_colors(s,:),...
                    sum(plot_onset_sbj{reg_ix}.(roi_list{roi_ix})==s),1);
            end
            violins(violin_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
            violins(violin_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
            violins(violin_ix).ScatterPlot.CData = scat_colors;
        else
            % Change violin color to match ROI
            violins(violin_ix).ViolinColor = roi_colors{roi_ix};
        end
    end
    
    % Add label and min RT for perspective
    ax = gca;
    ax.YLim = st.stat_lim;
    ax.YLabel.String = 'Time (s)';
    ax.Title.String = reg_names{reg_ix};
    set(gca,'FontSize',16);
    view([90 -90]);
    
    %% Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Error/results/ERP/GRP/onsets_ROI/'...
            model_id '/' stat_id '/' an_id '/' plt_id '/'];
        if ~exist(fig_dir,'dir')
            [~] = mkdir(fig_dir);
        end
        
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

end
