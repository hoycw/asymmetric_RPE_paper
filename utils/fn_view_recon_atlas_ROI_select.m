function fn_view_recon_atlas_ROI_select(SBJ, proc_id, plot_elecs, show_label,...
                                 hemi, atlas_id, roi_id, plot_roi, mirror, varargin)
error('didnt finish this! see run02_elec-Pipieline');
                             %% Plot a SBJ reconstruction with a subset of electrodes
% INPUTS:
%   SBJ [str] - ID of subject
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - gROI grouping to pick mesh and color specific ROIs
%       'LPFC','MPFC','OFC','INS','TMP','PAR'
%   plot_roi [str] - which surface mesh to plot
%   mirror [0/1] - plot the other hemi, 

%% Handle variables
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
        elseif strcmp(varargin{v},'fig_ftype')
            fig_ftype = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Define default options
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
if ~exist('fig_ftype','var')
    fig_ftype = 'fig';
end
if ~exist('save_fig','var')
    save_fig = 0;
end
if ~exist('mesh_alpha','var')
    % assume SEEG
    mesh_alpha = 0.3;
end

if show_label
    lab_arg = 'label';
else
    lab_arg = 'off';
end

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];    % MNI space
else
    reg_suffix = '';                % Patient space
end

if any(strcmp(plot_roi,{'deep','lat'}))
    [plot_roi_list, ~, ~] = fn_roi_label_styles(plot_roi);
else
    plot_roi_list = {plot_roi};
end
[~, ~, roi_field] = fn_roi_label_styles(roi_id);

[root_dir, ~] = fn_get_root_dir();
out_dir = [root_dir 'PRJ_Error/results/recons/'];

%% Load elec struct
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);
    cfgs = []; cfgs.channel = plot_elecs;
    elec = fn_select_elec(cfgs,elec);
    elec.color = fn_roi2color(elec.(roi_id));
    
    % Remove electrodes that aren't in atlas ROIs & hemisphere
    if mirror
        roi_elecs = fn_select_elec_lab_match(elec, 'b', atlas_id, roi_id);
        hemi_str = [hemi 'b'];
    else
        roi_elecs = fn_select_elec_lab_match(elec, hemi, atlas_id, roi_id);
        hemi_str = hemi;
    end
    good_elecs = intersect(roi_elecs, elec.label(any(plot_elecs,2)));
    % fn_select_elec messes up if you try to toss all elecs
    if isempty(good_elecs)
        elec = {};
        good_sbj(sbj_ix) = false;
    else
        cfgs = [];
        cfgs.channel = good_elecs;
        elec = fn_select_elec(cfgs, elec);
        
        % Mirror hemispheres
        if mirror
            elec.chanpos(~strcmp(elec.hemi,hemi),1) = ...
                -elec.chanpos(~strcmp(elec.hemi,hemi),1);
        end
        
        all_roi_labels = [all_roi_labels; elec.(roi_field)];
        all_roi_colors = [all_roi_colors; elec.color];
    end
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = ft_appendsens([],elec{good_sbj});
elec.roi   = all_roi_labels;    % appendsens strips that field
elec.color = all_roi_colors;    % appendsens strips that field

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
fig_name = [SBJ_id '_' atlas_id '_' roi_id '_' plot_roi '_' hemi_str '_' view_str];
h = figure('Name',fig_name);

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
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere',...
                 'facecolor', elec_tmp.color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(h, 'windowkeypressfcn',   @cb_keyboard);

%% Save figure
if save_fig
    saveas(gcf, [out_dir fig_name '.' fig_ftype]);
end

end
