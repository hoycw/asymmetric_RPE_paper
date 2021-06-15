function [roi_mesh,r_mesh,mtl_mesh] = fn_load_recon_mesh_ROI(atlas_id,plot_roi,hemi)
%% fn_load_recon_mesh_ROI creates pial surface meshes covering certain ROIs
% INPUTS:
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   plot_roi [str] - which surface mesh to plot
%       main options: 'lat','MPFC','deep','MTL','OFC'
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
% OUTPUTS:
%   roi_mesh [] - mesh surface of main ROI
%   r_mesh [] - if OFC, the second mesh for right hemisphere
%   mtl_mesh [] - if 'deep', mesh for the MTL (roi_mesh is insula)

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
roi_mesh = ft_prepare_mesh(cfg, seg);
if exist('r_seg','var')
    r_mesh = ft_prepare_mesh(cfg, r_seg);
else
    r_mesh = [];
end
if exist('mtl_seg','var')
    mtl_mesh = ft_prepare_mesh(cfg, mtl_seg);
else
    mtl_mesh = [];
end

end
