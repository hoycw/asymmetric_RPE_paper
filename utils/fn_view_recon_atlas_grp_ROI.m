function fn_view_recon_atlas_grp_ROI(SBJ_id, proc_id, atlas_id, roi_id, rcn, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - gROI grouping to pick mesh and color specific ROIs
%       'LPFC','MPFC','OFC','INS','TMP','PAR'
%   rcn [struct] - plotting options for the recon (i.e., recon_vars)
%       .plot_roi [str] - which surface mesh to plot
%       .reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%       .show_lab [0/1] - plot the electrode labels
%       .hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%       .mirror [0/1] - mirror elecs from one hemisphere to the other

%% Handle variables
if ~isstruct(rcn); error('rcn is not a struct!'); end
[root_dir, ~] = fn_get_root_dir();
out_dir = [root_dir 'PRJ_Error/results/recons/'];

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            rcn.view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            rcn.mesh_alpha = varargin{v+1};
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
if ~exist('save_fig','var');    save_fig = 0; end
if ~exist('fig_ftype','var');   fig_ftype = 'fig'; end

% ROI info
rcn = fn_process_recon_vars(rcn);

%% Load elec struct
SBJs = fn_load_SBJ_list(SBJ_id);

[elec_sbj, good_sbj] = fn_load_grp_elec_ROI(SBJs,proc_id,atlas_id,roi_id,rcn);

%% Combine elec structs
% Grab ROIs and colors to add back
all_roi_labels = {};
all_roi_colors = [];
for sbj_ix = 1:numel(SBJs)
    all_roi_labels = [all_roi_labels; elec_sbj{sbj_ix}.(rcn.roi_field)];
    all_roi_colors = [all_roi_colors; elec_sbj{sbj_ix}.color];
end
elec = ft_appendsens([],elec_sbj{good_sbj});
elec.roi   = all_roi_labels;    % appendsens strips that field
elec.color = all_roi_colors;    % appendsens strips that field

%% Load Surface Mesh based on ROI
[roi_mesh, roi_mesh_lab] = fn_load_recon_mesh_ROI(atlas_id,rcn.plot_roi,rcn.hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
fig_name = [SBJ_id '_' atlas_id '_' roi_id '_' rcn.plot_roi '_' rcn.hemi_str '_' rcn.view_str];
fig = fn_plot_recon_mesh(elec, roi_mesh, roi_mesh_lab, rcn, fig_name);

%% Save figure
if save_fig
    saveas(gcf, [out_dir fig_name '.' fig_ftype]);
end

end
