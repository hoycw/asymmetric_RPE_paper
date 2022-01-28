function fn_view_recon_atlas_grp_stat_ROI_LME_cat(SBJ_id, proc_id, an_id, model_id, stat_id, cat_id, reg_type, show_labels,...
                            hemi, atlas_id, roi_id, plot_roi, mirror, varargin)
%% Plot a reconstruction with electrodes colored according to statistics
%   This version uses categories based on EpnRPE model
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pcon': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%       skip_reg [str] - name of one regressor to zero out and skip (not plot)

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Process plotting params
% Error cases
if strcmp(hemi,'b') && ~strcmp(plot_roi,'OFC')
    error('hemi must be l or r for all non-OFC plots');
end
if ~any(strcmp(plot_roi,{'LPFC','MPFC','INS','OFC','TMP','PAR','lat','deep'}))
    error('roi_id needs to be a lobe, "lat", or "deep"');
end

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
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
% Add default view_angle if not defined
if ~exist('view_angle','var')
    view_angle = fn_get_view_angle(hemi,plot_roi);
    view_str = 'def';
end
% Adjust view angle if custom
if ischar(view_angle)
    view_str = view_angle;
    if (strcmp(hemi,'l') && strcmp(view_angle,'med')) || (strcmp(hemi,'r') && strcmp(view_angle,'lat'))
        view_angle = [90 0];
    elseif (strcmp(hemi,'l') && strcmp(view_angle,'lat')) || (strcmp(hemi,'r') && strcmp(view_angle,'med'))
        view_angle = [-90 0];
    end
end
if ~exist('mesh_alpha','var')
    mesh_alpha = 0.25;   % assume SEEG
end
if mirror && strcmp(hemi,'b')
    error('why mirror if hemi b?');
end
if ~exist('save_fig','var')
    save_fig = 0;
end
if ~exist('fig_ftype','var')
    fig_ftype = 'fig';
end
if show_labels
    lab_arg = 'label';
else
    lab_arg = 'off';
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

% ROI info
if any(strcmp(plot_roi,{'deep','lat'}))
    [plot_roi_list, ~] = fn_roi_label_styles(plot_roi);
else
    plot_roi_list = {plot_roi};
end
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
[cat_lab, cat_names, cat_colors, ~, ~] = fn_pnRPE_category_label_styles(cat_id);
venn_colors = fn_venn_colors(0,'model_lab',mdl.model_lab);
all_color   = [0.1 0.1 0.1];

SBJs = fn_load_SBJ_list(SBJ_id);

stats_fname = [root_dir 'PRJ_Error/data/GRP/stats/' model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname);

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

if strcmp(plot_roi,'deep')
    roi_label = 'INS';
elseif strcmp(plot_roi,'MPFC')
    roi_label = 'MPFC';
else
    error('SBJ08g not run yet for anything except MPFCINS');
end

%% Get significant channel labels
beta_roi_ix = strcmp(beta_chan.label,roi_label);
% category order from David's SBJ08g for [pRPE, nRPE, sRPE, uRPE]:
david_cat_order = {'pRPE','nRPE','uRPE','sRPE'};
sig_cat_elecs = cell(size(cat_lab));
for cat_ix = 1:length(cat_lab)
    david_ix = strcmp(david_cat_order,cat_lab{cat_ix});
    elec_ix = beta_chan.chancat_ix{beta_roi_ix}{david_ix,1};
    sig_cat_elecs{cat_ix} = beta_chan.chan_label{beta_roi_ix}(elec_ix);
    
    % strip extra SBJ label (!!! remove once this is fixed in SBJ08g!!!)
    for ch_ix = 1:length(sig_cat_elecs{cat_ix})
        sig_cat_elecs{cat_ix}{ch_ix} = sig_cat_elecs{cat_ix}{ch_ix}(6:end);
    end
end

%% Check ROI and significance matches
elec_sbj    = cell([numel(SBJs) 1]);
elec_sig    = cell([numel(SBJs) 1]);
good_sbj    = true([numel(SBJs) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    
    % Append SBJ name to labels
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
    end
    
    % Select ROI mesh matches
    plot_elecs = zeros([numel(elec_sbj{sbj_ix}.label) numel(plot_roi_list)]);
    for roi_ix = 1:numel(plot_roi_list)
        plot_elecs(:,roi_ix) = strcmp(elec_sbj{sbj_ix}.(roi_field),plot_roi_list{roi_ix});
    end
    % Remove electrodes that aren't in atlas ROIs & hemisphere
    if mirror
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, 'b', atlas_id, roi_id);
        hemi_str = [hemi 'b'];
    else
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
        hemi_str = hemi;
    end
    plot_roi_elecs = intersect(roi_elecs, elec_sbj{sbj_ix}.label(any(plot_elecs,2)));
    
    if ~isempty(plot_roi_elecs)
        % Mirror hemispheres
        if mirror
            elec_sbj{sbj_ix}.chanpos(~strcmp(elec_sbj{sbj_ix}.hemi,hemi),1) = ...
                -elec_sbj{sbj_ix}.chanpos(~strcmp(elec_sbj{sbj_ix}.hemi,hemi),1);
        end
        
        % Find significant elecs in ROI
        sig_cat_roi_elecs = intersect(plot_roi_elecs,vertcat(sig_cat_elecs{:}));
        if ~isempty(sig_cat_roi_elecs)
            cfgs = []; cfgs.channel = sig_cat_roi_elecs;
            elec_sig{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});
            fprintf('\t%s has %i sig channels in %s hemi %s\n',SBJ,size(sig_cat_roi_elecs,1),atlas_id,hemi);
        else
            % Print no significant elecs
            good_sbj(sbj_ix) = false;
            fprintf(2,'\t%s has %i channels in %s hemi %s, but none are significant\n',...
                SBJ,numel(plot_roi_elecs),atlas_id,hemi);
        end
    else
        % Print no ROI match
        good_sbj(sbj_ix) = false;
        fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJ,atlas_id,hemi_str);
    end
    
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = ft_appendsens([],elec_sig{good_sbj});
elec.color = nan([numel(elec.label) 3]);
for cat_ix = 1:numel(cat_lab)
    cat_idx = ismember(elec.label,sig_cat_elecs{cat_ix});
    elec.color(cat_idx,:) = repmat(cat_colors{cat_ix},sum(cat_idx),1);
end

%% Load Atlas
atlas = fn_load_recon_atlas([],atlas_id);

% Get Atlas-ROI mapping
atlas_labels = fn_atlas_roi_select_mesh(atlas_id, plot_roi, hemi);

% Can't plot unconnected meshes (I think), so create two meshes
if strcmp(plot_roi,'OFC')
    % Treat R hemi as new ROI
    r_ix = ~cellfun(@isempty,strfind(atlas_labels,'rh'));
    r_labels = atlas_labels(r_ix);
    atlas_labels = atlas_labels(~r_ix);
elseif strcmp(plot_roi,'deep')
    mtl_ix = ~cellfun(@isempty,strfind(atlas_labels,'Hippocampus')) | ...
             ~cellfun(@isempty,strfind(atlas_labels,'Amygdala'));
%     mtl_labels = atlas_labels(mtl_ix);
    atlas_labels = atlas_labels(~mtl_ix);
end

%% Select ROI mesh
cfg = [];
cfg.inputcoord = atlas.coordsys;
cfg.atlas = atlas;
cfg.roi = atlas_labels;
roi_mask = ft_volumelookup(cfg,atlas);
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = roi_mask;

if exist('r_labels','var')
    cfg.roi = r_labels;
    r_mask  = ft_volumelookup(cfg,atlas);
    
    r_seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    r_seg.brain = r_mask;
elseif exist('mtl_labels','var')
    cfg.roi  = mtl_labels;
    mtl_mask = ft_volumelookup(cfg,atlas);
    
    mtl_seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    mtl_seg.brain = mtl_mask;
end

cfg = [];
cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 100000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
roi_mesh = ft_prepare_mesh(cfg, seg);
if exist('r_seg','var')
    r_mesh = ft_prepare_mesh(cfg, r_seg);
elseif exist('mtl_seg','var')
    mtl_mesh = ft_prepare_mesh(cfg, mtl_seg);
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
if save_fig
    out_dir = [root_dir 'PRJ_Error/results/' an_dir '/GRP/recon_cat/' model_id '/' cat_id '/' stat_id '/' an_id '/'];
    if ~exist(out_dir,'dir')
        [~] = mkdir(out_dir);
    end
end
plot_name = [SBJ_id '_' model_id '_' cat_id '_' stat_id '_' an_id '_' plot_roi '_' hemi_str '_' view_str];
f = figure('Name',plot_name);

% Plot 3D mesh
ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
if exist('r_mesh','var')
    ft_plot_mesh(r_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
elseif exist('mtl_mesh','var')
    ft_plot_mesh(mtl_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
end

% Plot electrodes on top
for e = 1:numel(elec.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(f, 'windowkeypressfcn',   @cb_keyboard);

if save_fig
    fig_fname = [out_dir plot_name '.' fig_ftype];
    fig_fname = strrep(fig_fname,'*','x');
    saveas(f,fig_fname);
end

end
