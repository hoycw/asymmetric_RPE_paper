function SBJ08h_HFA_plot_grp_mLME_cat_recon(SBJ_id, proc_id, an_id, model_id, stat_id, cat_id,...
                            atlas_id, roi_id, rcn, varargin)
%% Plot a reconstruction with electrodes colored according to statistics
%   This version uses categories based on EpnRPE model
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   model_id [str] - ID of the model to load
%   stat_id [str] - ID of the stats (e.g., 'mLME_St0t6_WL05_WS25')
%   cat_id [str] - ID of the categories to group the coefficients
%       should be 'puns' for pos, neg, signed, unsigned/salience
%   atlas_id [str] - ID of atlas to select ROIs: {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%   rcn [struct] - options for recon plotting, see fn_process_recon_vars.m
%       .plot_roi [str] - which surface mesh to plot ('INS','MPFC',etc.)
%       .hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%       .mirror [0/1] - mirror elecs from one hemisphere to the other
%   OPTIONAL:
%       skip_reg [str] - name of one regressor to zero out and skip (not plot)

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Process plotting params
if ~isstruct(rcn); error('rcn is not a struct!'); end

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            rcn.view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            rcn.mesh_alpha = varargin{v+1};
        elseif strcmp(varargin{v},'skip_reg') && ischar(varargin{v+1})
            skip_reg = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype')
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Implement the default options
if ~exist('save_fig','var');    save_fig = 1; end
if ~exist('fig_ftype','var');   fig_ftype = 'fig'; end

% ROI info
rcn = fn_process_recon_vars(rcn);
[~, ~, roi_field] = fn_roi_label_styles(roi_id);

%% Load data
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
if strcmp(mdl.model_lab,'RL3D')
    [reg_lab, reg_names, reg_colors, reg_styles] = fn_performanceRL_regressor_label_styles(mdl.model_lab);
else
    [reg_lab, reg_names, reg_colors, reg_styles, reg_mrkrs] = fn_regressor_label_styles(mdl.model_lab);
end
if numel(reg_lab) < 2 || numel(reg_lab) > 3; error('why venn?'); end
[cat_lab, cat_names, cat_colors, ~, ~] = fn_puns_category_label_styles(cat_id);

SBJs = fn_load_SBJ_list(SBJ_id);

load([root_dir 'PRJ_Error/data/GRP/stats/' model_id '_' stat_id '_' an_id '_hfa_chancoef.mat']);

% Identify analysis type
if contains(an_id,'ERP')
    an_dir = 'ERP';
elseif contains(an_id,'HG')
    an_dir = 'HFA';
elseif  contains(an_id,'TFR')
    an_dir = 'TFR';
else
    error('unknown an_id');
end

if strcmp(rcn.plot_roi,'INS')
    roi_label = 'INS';
elseif strcmp(rcn.plot_roi,'MPFC')
    roi_label = 'MPFC';
else
    error('SBJ08g not run yet for anything except MPFCINS');
end

%% Get significant channel labels
beta_roi_ix = strcmp(beta_chan.label,roi_label);
if ~all(strcmp(beta_chan.chancat_label,cat_lab)); error('category label order mismatch'); end
sig_cat_elecs = cell(size(cat_lab));
for cat_ix = 1:length(cat_lab)
    elec_ix = beta_chan.chancat_ix{beta_roi_ix}{cat_ix,1};
    sig_cat_elecs{cat_ix} = beta_chan.chan_label{beta_roi_ix}(elec_ix);
end

%% Load elec and select within ROI and hemisphere
[elec_sbj, ~] = fn_load_grp_elec_ROI(SBJs,proc_id,atlas_id,roi_id,rcn);

%% Combine significant electrodes within the ROI and hemisphere
elec_sig = cell([numel(SBJs) 1]);
good_sbj = true([numel(SBJs) 1]);
for sbj_ix = 1:numel(SBJs)
    if ~isempty(elec_sbj{sbj_ix}.label)
        % Find significant elecs in ROI
        sig_cat_roi_elecs = intersect(elec_sbj{sbj_ix}.label,vertcat(sig_cat_elecs{:}));
        if ~isempty(sig_cat_roi_elecs)
            cfgs = []; cfgs.channel = sig_cat_roi_elecs;
            elec_sig{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});
            fprintf('\t%s has %i sig channels in %s (hemi %s)\n',SBJs{sbj_ix},size(sig_cat_roi_elecs,1),rcn.plot_roi,rcn.hemi);
        else
            % Print no significant elecs
            good_sbj(sbj_ix) = false;
            fprintf(2,'\t%s has %i channels in %s (hemi %s), but none are significant\n',...
                SBJs{sbj_ix},numel(elec_sbj{sbj_ix}.label),rcn.plot_roi,rcn.hemi);
        end
    else
        % Print no ROI match
        good_sbj(sbj_ix) = false;
%         fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJs{sbj_ix},atlas_id,rcn.hemi_str);
    end
    clear sig_cat_roi_elecs cfgs
end

%% Combine elec structs
elec = ft_appendsens([],elec_sig{good_sbj});
elec.color = nan([numel(elec.label) 3]);
for cat_ix = 1:numel(cat_lab)
    cat_idx = ismember(elec.label,sig_cat_elecs{cat_ix});
    elec.color(cat_idx,:) = repmat(cat_colors{cat_ix},sum(cat_idx),1);
end

%% Create 3D mesh based on atlas ROIs
[roi_mesh, roi_mesh_lab] = fn_load_recon_mesh_ROI(atlas_id,rcn.plot_roi,rcn.hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
if save_fig
    out_dir = [root_dir 'PRJ_Error/results/' an_dir '/GRP/recon_cat/' model_id '/' cat_id '/' stat_id '/' an_id '/'];
    if ~exist(out_dir,'dir'); [~] = mkdir(out_dir); end
end
plot_name = [SBJ_id '_' model_id '_' cat_id '_' stat_id '_' an_id '_' rcn.plot_roi '_' rcn.hemi_str '_' rcn.view_str];

fig = fn_plot_recon_mesh(elec, roi_mesh, roi_mesh_lab, rcn, plot_name);

if save_fig
    fig_fname = [out_dir plot_name '.' fig_ftype];
    fig_fname = strrep(fig_fname,'*','x');
    saveas(fig,fig_fname);
end

end
