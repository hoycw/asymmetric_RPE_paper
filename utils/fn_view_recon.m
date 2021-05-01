function fn_view_recon(SBJ, proc_id, plot_type, view_space, reg_type,...
                       show_labels, hemi, plot_out, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of processing pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Variable Handline
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('view_angle','var')
    view_angle     = [-90 0];
end
if ~exist('mesh_alpha','var')
    if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
        mesh_alpha = 0.3;
    else
        mesh_alpha = 0.8;
    end
end
if show_labels
    lab_arg = 'label';
else
    lab_arg = 'off';
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];    % MNI space
else
    reg_suffix = '';                % Patient space
end

%% Load elec struct
if isempty(proc_id)
    % Original elec files
    elec_fname = eval(['SBJ_vars.recon.elec_' view_space reg_suffix]);
    slash = strfind(elec_fname,'/'); elec_suffix = elec_fname(slash(end)+numel(SBJ)+2:end-4);
    
    tmp = load(elec_fname);
    elec_var_name = fieldnames(tmp);
    if ~strcmp(elec_var_name,elec_suffix)
        warning(['\t!!!! ' SBJ ' elec names in variable and file names do not match! file=' elec_suffix '; var=' elec_var_name{1}]);
    end
    eval(['elec = tmp.' elec_var_name{1} ';']); clear tmp;
else
    % Preprocessed (bipolar) elec files
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_',view_space,reg_suffix,'.mat']);
end

%% Remove electrodes that aren't in hemisphere
if ~plot_out
    cfgs = [];
    cfgs.channel = fn_select_elec_lab_match(elec, hemi, [], []);
    elec = fn_select_elec(cfgs, elec);
end

%% Load brain recon
if strcmp(plot_type,'3d')
    mesh = fn_load_recon_mesh(SBJ,view_space,reg_type,hemi);
elseif strcmp(plot_type,'ortho')
    mri = fn_load_recon_mri(SBJ,view_space,reg_type);
else
    error(['Unknown plot_type: ' plot_type]);
end

%% Orthoplot (pat/mni, v only, 0/1 labels)
if strcmp(plot_type,'ortho')
    % ft_electrodeplacement only plots elec.elecpos, so swap in chanpos
    elec.elecpos = elec.chanpos;
    if isfield(elec,'tra')
        elec = rmfield(elec, 'tra');    % error if elec.tra shows the difference between original elecpos and new chanpos post-reref
    end
    cfg = [];
    cfg.elec = elec;
    ft_electrodeplacement(cfg, mri);
end

%% 3D Surface + Grids (3d, pat/mni, v/s, 0/1)
if strcmp(plot_type,'3d')
    h = figure;
    
    % Plot 3D mesh
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', lab_arg);
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(h, 'windowkeypressfcn',   @cb_keyboard);
end

end
