%% Run scripts to get data for Thesis Table 1
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);

proc_id = 'main_ft';

%% Load Paradigm Variables
% prdm_vars = cell(size(SBJs));
% fnames    = cell(size(SBJs));
% itis      = cell(size(SBJs));
% blk_order = cell(size(SBJs));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('%s:\n',SBJ);
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    
    fprintf('\t%s srates: ',SBJs{s}); fprintf('%.2f, ',itis{s}); fprintf('\n');
    clear SBJ SBJ_vars
end


