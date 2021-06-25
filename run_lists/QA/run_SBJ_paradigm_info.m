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
prdm_vars = cell(size(SBJs));
fnames    = cell(size(SBJs));
itis      = cell(size(SBJs));
blk_order = cell(size(SBJs));
for s = 1:numel(SBJs)
    prdm_vars{s} = load([root_dir 'PRJ_Error/data/' SBJs{s} '/03_events/' SBJs{s} '_prdm_vars.mat']);
    fnames{s} = fieldnames(prdm_vars{s});
    load([root_dir 'PRJ_Error/data/' SBJs{s} '/03_events/' SBJs{s} '_bhv_' proc_id '_final.mat']);
    new_blk_ix = [1; find(diff(bhv.blk)~=0)+1];
    blk_order{s} = bhv.cond(new_blk_ix);
    itis{s} = unique(bhv.ITI_type);
    
    clear bhv
end
all_fnames = unique(vertcat(fnames{:}));

%% Compare Paradigm Variables
% block order
for s = 1:numel(SBJs)
    fprintf('\t%s ITIs: ',SBJs{s}); fprintf('%.2f, ',itis{s}); fprintf('\n');
end

% block order
for s = 1:numel(SBJs)
    fprintf('%s block order: ',SBJs{s}); fprintf('%s, ',blk_order{s}{:}); fprintf('\n');
end

% tolerance limits
for s = 1:numel(SBJs)
    fprintf('%s tolerance limits: ',SBJs{s}); fprintf('%.3f, ',prdm_vars{s}.tol_lim); fprintf('\n');
end


