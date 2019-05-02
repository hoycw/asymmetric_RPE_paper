function atlas = fn_load_recon_atlas(SBJ,atlas_id)
%% Load an atlas for selecting ROIs from recon MRI

fprintf('Using atlas: %s\n',atlas_id);
[root_dir, ~] = fn_get_root_dir();

if ~isempty(SBJ)
    % Individual SBJ atlases
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    if strcmp(atlas_id,'DK')
        atlas      = ft_read_atlas(SBJ_vars.recon.fs_DK); % Desikan-Killiany (+volumetric)
        atlas.coordsys = 'acpc';
    elseif strcmp(atlas_id,'Dx')
        atlas      = ft_read_atlas(SBJ_vars.recon.fs_Dx); % Destrieux (+volumetric)
        atlas.coordsys = 'acpc';
    else
        error(['altas: ' atlas_id ' not compatible with specific SBJ']);
    end
else
    fscolin_dir = [root_dir 'PRJ_Error/data/atlases/freesurfer/fscolin/'];
    if strcmp(atlas_id,'DK')
        atlas      = ft_read_atlas([fscolin_dir 'fscolin_aparc+aseg.mgz']); % Desikan-Killiany (+volumetric)
        atlas.coordsys = 'acpc';
    elseif strcmp(atlas_id,'Dx')
        atlas      = ft_read_atlas([fscolin_dir 'fscolin_aparc.a2009s+aseg.mgz']); % Destrieux (+volumetric)
        atlas.coordsys = 'acpc';
    elseif strcmp(atlas_id,'Yeo7')
        atlas = fn_read_atlas(atlas_id);
        atlas.coordsys = 'mni';
    elseif strcmp(atlas_id,'Yeo17')
        atlas = fn_read_atlas(atlas_id);
        atlas.coordsys = 'mni';
    else
        error(['atlas_id unknown: ' atlas_id]);
    end    
end

atlas.name = atlas_id;

end
