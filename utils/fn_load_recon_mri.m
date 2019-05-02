function mri = fn_load_recon_mri(SBJ,view_space,reg_type)
%% Load the surface mesh of a recon
%   SBJ [str] - subject ID to plot
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(view_space,'pat')
    mri = ft_read_mri(SBJ_vars.recon.fs_T1);
    mri.coordsys = 'acpc';
elseif strcmp(view_space,'mni')
    if strcmp(reg_type,'v')
        mri = ft_read_mri([ft_dir 'template/anatomy/single_subj_T1_1mm.nii']);
        mri.coordsys = 'mni';
    elseif strcmp(reg_type,'s')
        error('ortho plot with surface based registration doesnt make sense!');
    end
else
    error(['Unknown view_space: ' view_space]);
end


end
