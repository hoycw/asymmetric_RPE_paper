function [roi_mesh, roi_mesh_lab] = fn_load_recon_mesh_ROI(atlas_id,plot_roi,hemi)
%% fn_load_recon_mesh_ROI creates pial surface meshes covering certain ROIs
% INPUTS:
%   atlas_id [str] - ID of atlas to select ROIs: {'DK','Dx','Yeo7','Yeo17'}
%   plot_roi [str] - which surface mesh to plot
%       {'LPFC','MPFC','INS','OFC','TMP','PAR','MTL','lat','deep'}
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
% OUTPUTS:
%   roi_mesh [cell array] - mesh surfaces in order: (main ROI, right ROI, MTL)
%       if OFC, the second mesh is for the right hemisphere
%       if 'deep', second mesh is for the MTL (main ROI is insula)
%   roi_mesh_lab [cell array] - names of the meshes
%       first should be 'ROI' (whatever is plot_roi)
%       second should be 'r_ROI' or 'MTL'

%% Load atlas and select subregion labels based on ROI to be plotted
atlas = fn_load_recon_atlas([],atlas_id);

% Get Atlas-ROI mapping
if strcmp(plot_roi,'MPFCINS')
    ins_atlas_labels = fn_atlas_roi_select_mesh(atlas_id, 'INS', hemi);
    if strcmp(hemi,'r'); other_hemi = 'l'; else other_hemi = 'r'; end
    mpfc_atlas_labels = fn_atlas_roi_select_mesh(atlas_id, 'MPFC', other_hemi);
    atlas_labels = [mpfc_atlas_labels; ins_atlas_labels];
else
    atlas_labels = fn_atlas_roi_select_mesh(atlas_id, plot_roi, hemi);
end

% Can't plot unconnected meshes (I think), so create two meshes
if any(strcmp(plot_roi,{'OFC','MPFCINS'}))
    % Treat R hemi as new ROI
    r_ix = contains(atlas_labels,'rh');
    r_labels = atlas_labels(r_ix);
    atlas_labels = atlas_labels(~r_ix);
elseif strcmp(plot_roi,'deep')
    mtl_ix = contains(atlas_labels,'Hippocampus') | contains(atlas_labels,'Amygdala');
    mtl_labels = atlas_labels(mtl_ix);
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
cfg.numvertices = 200000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
roi_mesh{1} = ft_prepare_mesh(cfg, seg);
roi_mesh_lab = {plot_roi};
if exist('r_seg','var')
    roi_mesh{length(roi_mesh)+1} = ft_prepare_mesh(cfg, r_seg);
    roi_mesh_lab{length(roi_mesh_lab)+1} = ['r_' plot_roi];
end
if exist('mtl_seg','var')
    roi_mesh{length(roi_mesh)+1} = ft_prepare_mesh(cfg, mtl_seg);
    roi_mesh_lab{length(roi_mesh_lab)+1} = 'MTL';
end

end