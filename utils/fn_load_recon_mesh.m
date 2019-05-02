function mesh = fn_load_recon_mesh(SBJ,view_space,reg_type,hemi)
%% Load the surface mesh of a recon
%   SBJ [str] - subject ID to plot
%       '' if view_space = 'mni'
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if ~isempty(SBJ)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
end

if strcmp(view_space,'pat')
    if strcmp(hemi,'r') || strcmp(hemi,'l')
        mesh = ft_read_headshape([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_' hemi 'h.mat']);
    elseif strcmp(hemi,'b')
        mesh = ft_read_headshape({[SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_rh.mat'],...
                                    [SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_lh.mat']});
    else
        error(['Unknown hemisphere selected: ' hemi]);
    end
    mesh.coordsys = 'acpc';
elseif strcmp(view_space,'mni')
    if strcmp(reg_type,'v')
        if strcmp(hemi,'r')
            load([ft_dir 'template/anatomy/surface_pial_right.mat']);
        elseif strcmp(hemi,'l')
            load([ft_dir 'template/anatomy/surface_pial_left.mat']);
        elseif strcmp(hemi,'b')
            load([ft_dir 'template/anatomy/surface_pial_both.mat']);
        else
            error(['Unknown hemisphere option: ' hemi]);
        end
%         mesh.coordsys = 'mni';
    elseif strcmp(reg_type,'s')
        if strcmp(hemi,'r') || strcmp(hemi,'l')
            mesh = ft_read_headshape([root_dir 'PRJ_Error/data/atlases/freesurfer/fsaverage/' hemi 'h.pial']);
        elseif strcmp(hemi,'b')
            error('hemisphere "b" not yet implemented for reg_type: "srf"!');
            mesh = ft_read_headshape([ft_dir 'subjects/fsaverage/surf/' hemi 'h.pial']);
        else
            error(['Unknown hemisphere option: ' hemi]);
        end
        mesh.coordsys = 'fsaverage';
    else
        error(['Unknown registration type (reg_type): ' reg_type]);
    end
else
    error(['Unknown view_space: ' view_space]);
end

end
