function SBJ06a_ERP_stats(SBJ,conditions,proc_id,an_id)
% Calculates ERPs, computes cluster-based statistics, and plots the results

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
[cond_lab, ~, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(conditions,bhv);

% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Compute ERPs
% Preprocess the data
cfgpp = [];
cfgpp.hpfilter  = hp_yn;
cfgpp.hpfreq    = hp_freq;
cfgpp.hpfiltord = hp_filtord;            % Leaving blank causes instability error, 1 or 2 works 
cfgpp.lpfilter  = lp_yn;
cfgpp.lpfreq    = lp_freq;
roi = ft_preprocessing(cfgpp,roi);

% Cut into Trials
if strcmp(an.event_type,'stim')
    events = bhv.trl_onset;
elseif strcmp(event_type,'resp')
    events = bhv.rsp_onset;
elseif strcmp(event_type,'fb')
    events = bhv.fb_onset;
else
    error(stract('ERROR: unknown event_type ',event_type));
end
roi_trl = fn_ft_cut_trials_equal_len(roi,events,cond_idx,trial_lim_s*roi.fsample);

% Baseline Correction
cfg = [];
cfg.demean         = demean_yn;
cfg.baselinewindow = bsln_lim;
roi_trl = ft_preprocessing(cfg,roi_trl);

% Average ERPs
roi_erp = {};
n_trials = zeros([1 numel(cond_lab)]);
cfgavg = [];
cfgavg.keeptrials = 'yes';
for cond_ix = 1:numel(cond_lab)
    cfgavg.trials = find(cond_idx==cond_ix);
    roi_erp{cond_ix} = ft_timelockanalysis(cfgavg,roi_trl);
    % Grab n_trials for design matrix
    n_trials(cond_ix) = size(roi_erp{cond_ix}.trial,1);
end

%% Run Statistics
% Create design matrix
% design = zeros(2,size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1));
% % Conditions (Independent Variable)
% design(1,1:size(roi_erp_con.trial,1)) = 1;
% design(1,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
% % Trial Numbers
% design(2,1:size(roi_erp_con.trial,1)) = 1;
% design(2,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
design = zeros(2,sum(n_trials));
for cond_ix = 1:numel(cond_lab)
    if cond_ix==1
        design(1,1:n_trials(cond_ix)) = cond_ix;                                % Conditions (Independent Variable)
        design(2,1:n_trials(cond_ix)) = 1:n_trials(cond_ix);                    % Trial Numbers
    else
        design(1,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= cond_ix; % Conditions (Independent Variable)
        design(2,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= 1:n_trials(cond_ix);
    end
end

% Prepare neighbors layout
% cfgn = [];
% cfgn.method  = 'distance';
% cfgn.layout  = 'ordered';
% cfgn.channel = elecs;
% neighbors    = ft_prepare_neighbours(cfgn,roi_erp_allch{1});
% for ch_ix = 1:numel(roi_erp{1}.label)
%     neighbors(ch_ix).label = roi_erp{1}.label{ch_ix};
%     neighbors(ch_ix).neighblabel = {};
% end

% Calculate statistics
cfg_stat.design           = design;
[stat] = ft_timelockstatistics(cfg_stat, roi_erp{:});

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_ROI_',conditions,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','roi_erp','stat');

end
