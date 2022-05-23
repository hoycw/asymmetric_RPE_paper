function fn_view_recon(SBJ, proc_id, rcn)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of processing pipeline, used to pick elec file
%   rcn [struct] - options for recon plotting, see fn_process_recon_vars.m

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Variable Handline
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

rcn = fn_process_recon_vars(rcn);
if ~any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
    rcn.mesh_alpha = 0.8;
end

%% Load elec struct
if isempty(proc_id)
    % Original elec files
    elec_fname = eval(['SBJ_vars.recon.elec_' rcn.view_space rcn.reg_suffix]);
    slash = strfind(elec_fname,'/'); elec_suffix = elec_fname(slash(end)+numel(SBJ)+2:end-4);
    
    tmp = load(elec_fname);
    elec_var_name = fieldnames(tmp);
    if ~strcmp(elec_var_name,elec_suffix)
        warning(['\t!!!! ' SBJ ' elec names in variable and file names do not match! file=' elec_suffix '; var=' elec_var_name{1}]);
    end
    eval(['elec = tmp.' elec_var_name{1} ';']); clear tmp;
else
    % Preprocessed (bipolar) elec files
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_',rcn.view_space,rcn.reg_suffix,'.mat']);
end

%% Remove electrodes that aren't in hemisphere
if ~rcn.plot_out
    cfgs = [];
    cfgs.channel = fn_select_elec_lab_match(elec, rcn.hemi, [], []);
    elec = fn_select_elec(cfgs, elec);
end

%% Load brain recon
if strcmp(rcn.plot_type,'3d')
    brain_mesh = fn_load_recon_mesh(SBJ,rcn.view_space,rcn.reg_type,'pial',rcn.hemi);
elseif strcmp(rcn.plot_type,'ortho')
    mri = fn_load_recon_mri(SBJ,rcn.view_space,rcn.reg_type);
else
    error(['Unknown plot_type: ' rcn.plot_type]);
end

%% Orthoplot (pat/mni, v only, 0/1 labels)
if strcmp(rcn.plot_type,'ortho')
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
if strcmp(rcn.plot_type,'3d')
    fig = fn_plot_recon_mesh(elec, brain_mesh, 'brain', rcn, '');
end

end
