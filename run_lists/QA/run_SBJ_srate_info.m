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
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('%s:\n',SBJ);
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    for r_ix = 1:numel(SBJ_vars.raw_file)
        % Load the data
        if strcmp(SBJ_vars.raw_file{r_ix}(end-2:end),'mat')
            load(SBJ_vars.dirs.raw_filename{r_ix});
        else
            cfg = [];
            cfg.dataset = SBJ_vars.dirs.raw_filename{r_ix};
            cfg.continuous = 'yes';
            cfg.channel = 'all';
            data = ft_preprocessing(cfg);
        end
        
        fprintf('\tR%d srate = %.1f\n',r_ix,data.fsample);
    end
    
    clear SBJ SBJ_vars data
end

%% Event Channel sampling rates for single unit patients
evnt_srates  = nan(size(SBJs));
macro_srates = nan(size(SBJs));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    if isfield(SBJ_vars.dirs,'nlx')
        evnt = fn_read_neuralynx_interp({[SBJ_vars.dirs.nlx{1} 'photo/' ...
            SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.nlx_suffix{1} '.ncs']});
        
        macro_fname = [SBJ_vars.dirs.nlx{1} 'macro/' SBJ_vars.ch_lab.nlx_nk_align{1}...
            SBJ_vars.ch_lab.nlx_suffix{1} '.ncs'];
        macro = fn_read_neuralynx_interp({macro_fname});
        
        evnt_srates(s) = evnt.fsample;
        macro_srates(s) = macro.fsample;
    end
    
    clear SBJ SBJ_vars evnt macro
end

% Print srates
for s = 1:numel(SBJs)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    if isfield(SBJ_vars.dirs,'nlx')
        fprintf('%s:\n',SBJs{s});
        fprintf('\tphoto srate = %.1f\n',evnt_srates(s));
        fprintf('\tmacro srate = %.1f\n',macro_srates(s));
    end
    clear SBJ_vars
end