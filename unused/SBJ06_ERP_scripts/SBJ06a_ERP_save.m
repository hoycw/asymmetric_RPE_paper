function SBJ06a_ERP_save(SBJ,proc_id,an_id)
%% Compute single trial ERPs
%   Select ROI channels, filter, cut to trials, baseline demean

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%%
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_bhv_',proc_id,'_final.mat'));

% Select Conditions of Interest
[cond_lab, ~, ~] = fn_condition_label_styles('DifFB'); % keep all conditions
cond_idx = fn_condition_index(cond_lab,bhv);

% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Preprocess the data
cfgpp = [];
cfgpp.hpfilter  = an.hp_yn;
cfgpp.hpfreq    = an.hp_freq;
cfgpp.hpfiltord = an.hp_filtord;            % Leaving blank causes instability error, 1 or 2 works 
cfgpp.lpfilter  = an.lp_yn;
cfgpp.lpfreq    = an.lp_freq;
roi = ft_preprocessing(cfgpp,roi);

%% Cut into Trials
%   for now, only baseline to the same time lock event (e.g., no stim
%   baseline for feedback-locked ERPs)
% Check that baseline will be included in data cut to trial_lim_s
if an.trial_lim_s(1) < an.bsln_lim(1)
    error(['ERROR: an.trial_lim_s does not include an.bsln_lim for an_id = ' an_id]);
end
if strcmp(an.evnt_lab,'S')
    events = bhv.trl_onset;
elseif strcmp(an.evnt_lab,'R')
    events = bhv.rsp_onset;
elseif strcmp(an.evnt_lab,'F')
    events = bhv.fb_onset;
else
    error(stract('ERROR: unknown an.evnt_lab ',an.evnt_lab));
end
roi_trl = fn_ft_cut_trials_equal_len(roi,events,cond_idx,an.trial_lim_s*roi.fsample);

%% Baseline Correction
cfg = [];
cfg.demean         = an.demean_yn;
cfg.baselinewindow = an.bsln_lim;
erp_trl = ft_preprocessing(cfg,roi_trl);

%% Save Results
data_out_fname = [SBJ_vars.dirs.proc SBJ '_ROI_' proc_id '_' an_id '.mat'];
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','erp_trl');

end
