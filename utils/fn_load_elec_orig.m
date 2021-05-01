function [elec] = fn_load_elec_orig(SBJ,view_space,reg_type)
%% function [elec] = fn_load_elec_orig()
% INPUTS:
%   SBJ [str]
%   view_space [str] - {'pat','mni'}
%   reg_type [str] - {'v','s',''};

[root_dir, ~] = fn_get_root_dir();
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];    % MNI space
else
    reg_suffix = '';                % Patient space
end

% Get file name
elec_fname = eval(['SBJ_vars.recon.elec_' view_space reg_suffix]);
slash = strfind(elec_fname,'/'); elec_suffix = elec_fname(slash(end)+numel(SBJ)+2:end-4);

% Load file
tmp = load(elec_fname);

% Rename file
elec_var_name = fieldnames(tmp);
if ~strcmp(elec_var_name,elec_suffix)
    warning(['\t!!!! ' SBJ ' elec names in variable and file names do not match! file=' elec_suffix '; var=' elec_var_name{1}]);
end
eval(['elec = tmp.' elec_var_name{1} ';']); clear tmp;

end
