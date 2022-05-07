function rcn = fn_process_recon_vars(rcn)
%% Process options for plotting recons
% INPUT Fields:
%  -Necessary:
%   .hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%  -Optional: (* indicates default)
%   .view_space [str] - 'pat' for patient, *'mni' for group
%   .plot_type [str] - 'ortho' for 3 orthogonal planes, *'3d' for surface meshes
%   .plot_out [0/1] - plot electrodes outside the brain, hemisphere, or ROI
%   .reg_type [str] - {'', 'v', 's'} choose volume-based or surface-based registration
%       default = 'v' unless view_space = 'pat', in which case default is ''
%   .show_lab [0/1] - plot the electrode labels (default = 0)
%   .plot_roi [str] - which surface mesh to plot (default = '')
%   .mirror [0/1] - mirror elecs from one hemisphere to the other
%       default = 1 (yes) unless plot_roi = 'OFC'
% OUTPUT Fields:
%   .view_angle [1x2 float] - angel to view the recon
%   .view_str [str] - label to identify the view when saving the figure
%   .hemi_str [str] - label of hemisphere (including mirroring) for saving
%   .mesh_alpha [float] - transparency from 0 (clear) to 1 (opaque)
%       deafult = 0.25 for SEEG
%   .lab_arg [str] - transformed 0/1 to 'off'/'label' for ft_plot_sens
%   .reg_suffix [str] - reg_type transformed for saving output
%   .plot_roi_list [cell array] - list of roi_lab that fall within plot_roi

%% Fill missing defaults
if ~isfield(rcn,'hemi'); error('rcn.hemi is the one required field!'); end
if ~isfield(rcn,'view_space');  rcn.view_space = 'mni'; end
if ~isfield(rcn,'plot_type');   rcn.plot_type = '3d'; end
if ~isfield(rcn,'plot_out');    rcn.plot_out = 0; end
if ~isfield(rcn,'show_lab');    rcn.show_lab = 0; end
if ~isfield(rcn,'plot_roi');    rcn.plot_roi = ''; end
if ~isfield(rcn,'mesh_alpha');  rcn.mesh_alpha = 0.25; end   % assume SEEG
if ~isfield(rcn,'reg_type')
    if strcmp(rcn.view_space,'pat'); rcn.reg_type = ''; else; rcn.reg_type = 'v'; end
end
if ~isfield(rcn,'mirror')
    if any(strcmp(rcn.plot_roi,{'OFC',''})); rcn.mirror = 0; else; rcn.mirror = 1; end
end

%% Error cases
if rcn.mirror && strcmp(rcn.hemi,'b'); error('why mirror if hemi b?'); end
if strcmp(rcn.hemi,'b') && ~any(strcmp(rcn.plot_roi,{'OFC',''}))
    error('hemi must be l or r for all non-OFC plots');
end
if ~any(strcmp(rcn.plot_roi,{'LPFC','MPFC','INS','OFC','TMP','PAR','MTL','lat','deep',''}))
    error('roi_id needs to be a lobe, "MTL", "lat", or "deep"');
end

%% View Angle
% Add default view_angle if not defined
if ~isfield(rcn,'view_angle')
    rcn.view_angle = fn_get_view_angle(rcn.hemi,rcn.plot_roi);
    rcn.view_str = 'def';
end
% Adjust view angle if custom
if ischar(rcn.view_angle)
    rcn.view_str = rcn.view_angle;
    if (strcmp(rcn.hemi,'l') && strcmp(rcn.view_angle,'med')) || (strcmp(rcn.hemi,'r') && strcmp(rcn.view_angle,'lat'))
        rcn.view_angle = [90 0];
    elseif (strcmp(rcn.hemi,'l') && strcmp(rcn.view_angle,'lat')) || (strcmp(rcn.hemi,'r') && strcmp(rcn.view_angle,'med'))
        rcn.view_angle = [-90 0];
    end
end

%% Convert plot_roi into list of general or specific ROIs
if any(strcmp(rcn.plot_roi,{'deep','lat'}))
    [rcn.plot_roi_list, ~] = fn_roi_label_styles(rcn.plot_roi);
else
    rcn.plot_roi_list = {rcn.plot_roi};
end

%% String handling
if rcn.show_lab
    rcn.lab_arg = 'label';
else
    rcn.lab_arg = 'off';
end
if strcmp(rcn.reg_type,'v') || strcmp(rcn.reg_type,'s')
    rcn.reg_suffix = ['_' rcn.reg_type];
else
    rcn.reg_suffix = '';
end
if rcn.mirror; rcn.hemi_str = [rcn.hemi 'b']; else; rcn.hemi_str = rcn.hemi; end

end