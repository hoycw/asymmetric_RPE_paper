function trl_info_cln = SBJ04_compile_clean_behavior(SBJ, proc_id, save_it)
% Reject all bad trials based on behavior and visual cleaning, then combine
%   data from multiple runs
% Criteria:
%   1. Bob bad epoch
%   2. Bad trial (interruption, corrupt data, etc.)
%   3. Bad behavior (no response, RT outlier, etc.)
%   4. Training trials
% Inputs:
%   SBJ [str]- the dataset to process (e.g., 'IR54')
%   proc_id [str] - name of processing pipeline (proc_vars)
%       should contain RT_std_thresh [int]- # standard deviation from mean for RT to be tossed as outlier
%   save_it [0/1] - save the cleaned up trl_info
% Outputs:
%   trl_info_clean [struct]- saves out final version after tossing all bad trials

% Directories
[root_dir, ~] = fn_get_root_dir();
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% Load data
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_proc_vars.m']);
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_bad_epochs_preproc.mat'));

% Get sampling rate
if SBJ_vars.low_srate(b_ix)~=0
    data_srate = SBJ_vars.low_srate(b_ix);
else
    data_srate = proc_vars.resample_freq;
end

% Load different trl_infos
trl_infos = cell(size(SBJ_vars.block_name));
trl_cnt   = zeros(size(SBJ_vars.block_name));
blk_cnt   = zeros(size(SBJ_vars.block_name));
blk_lens  = zeros(size(SBJ_vars.block_name));
blk_times = zeros(size(SBJ_vars.block_name));
for b_ix = 1:numel(SBJ_vars.block_name)
    % Get block timing
    if numel(SBJ_vars.raw_file)>1
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    end
    
    % Load trl_info
    tmp = load([SBJ_vars.dirs.events SBJ '_trl_info_auto' block_suffix '.mat']);
    trl_infos{b_ix} = tmp.trl_info;
    
    % Get block length
    tmp = load(strcat(SBJ_vars.dirs.import,SBJ,'_',...
        num2str(data_srate),'hz',block_suffix,'.mat'));
    blk_lens(b_ix) = size(tmp.data.trial{1},2);
    blk_times(b_ix) = tmp.data.time{1}(end);
    
    % Add run stats (trial, block, run counts)
    trl_cnt(b_ix) = numel(trl_infos{b_ix}.trl_n);
    blk_cnt(b_ix) = numel(setdiff(unique(trl_infos{b_ix}.blk),0));
    trl_infos{b_ix}.run = repmat(b_ix,size(trl_infos{b_ix}.trl_n));
    
    % Check if same fields
    if b_ix==1
        ti_fields = fieldnames(trl_infos{b_ix});
        prdm_fields = fieldnames(trl_infos{b_ix}.prdm);
    else
        if numel(fieldnames(trl_infos{b_ix}))~=numel(ti_fields) || ~all(strcmp(fieldnames(trl_infos{b_ix}),ti_fields))
            error(['Mismatched trl_info fields between ' SBJ_vars.block_name{b_ix} ' and ' SBJ_vars.block_name{1}]);
        end
        if numel(fieldnames(trl_infos{b_ix}.prdm))~=numel(prdm_fields) || ~all(strcmp(fieldnames(trl_infos{b_ix}.prdm),prdm_fields))
            error(['Mismatched prdm fields between ' SBJ_vars.block_name{b_ix} ' and ' SBJ_vars.block_name{1}]);
        end
        for p_ix = 1:numel(prdm_fields)
            if ischar(trl_infos{b_ix}.prdm.(prdm_fields{p_ix}))
                if ~strcmp(trl_infos{b_ix}.prdm.(prdm_fields{p_ix}),trl_infos{1}.prdm.(prdm_fields{p_ix}))
                    error(['Mismatch prdm field: ' prdm_fields{p_ix}]);
                end
            elseif ~all(trl_infos{b_ix}.prdm.(prdm_fields{p_ix})==trl_infos{1}.prdm.(prdm_fields{p_ix}))
                error(['Mismatch prdm field: ' prdm_fields{p_ix}]);
            end
        end
    end
end

%% Concatenate all fields
trl_info = trl_infos{1};
trl_info.ignore_trials = {trl_info.ignore_trials};
% Keep the original trl_info structs from each run
if numel(trl_infos)>1; trl_info.run_trl_info = trl_infos; end  

% Add properties of the individual blocks
trl_info.run_name = SBJ_vars.block_name;
trl_info.run_len  = blk_lens;
trl_info.run_time = blk_times;
% Combine trl_info structs if necessary
for b_ix = 2:numel(SBJ_vars.block_name)
    for f_ix = 1:numel(ti_fields)
        if numel(trl_infos{b_ix}.(ti_fields{f_ix}))==trl_cnt(b_ix)
            % Concatenate fields that don't need modification
            %   NOTE: keeping blk in original numbers (helps toss training)
            if any(strcmp(ti_fields{f_ix},{'cond','blk','run','hit','rt','tol','blk_trl_n','score','ITI','ITI_type'}))
                trl_info.(ti_fields{f_ix}) = vertcat(trl_info.(ti_fields{f_ix}),trl_infos{b_ix}.(ti_fields{f_ix}));
            % Modify then concatenate counts and indices
            elseif any(strcmp(ti_fields{f_ix},{'trl_onset','rsp_onset','fb_onset'}))
                trl_info.(ti_fields{f_ix}) = vertcat(trl_info.(ti_fields{f_ix}),trl_infos{b_ix}.(ti_fields{f_ix})+sum(blk_lens(1:b_ix-1)));
            elseif strcmp(ti_fields{f_ix},'trl_n')
                trl_info.trl_n = vertcat(trl_info.trl_n,trl_infos{b_ix}.trl_n+sum(trl_cnt(1:b_ix-1)));
            elseif strcmp(ti_fields{f_ix},'time')
                trl_info.time = vertcat(trl_info.time,trl_infos{b_ix}.time+sum(blk_times(1:b_ix-1)));
            else
                error(['Unknown field: ' ti_fields{f_ix}]);
            end
        elseif strcmp(ti_fields{f_ix},'ignore_trials')
            trl_info.ignore_trials = [trl_info.ignore_trials {trl_infos{b_ix}.ignore_trials}];
        elseif any(strcmp(ti_fields{f_ix},{'SBJ','rt_type'}))
            if trl_info.(ti_fields{f_ix})~=trl_infos{b_ix}.(ti_fields{f_ix})
                error([ti_fields{f_ix} ' mismatch!']);
            end
        elseif strcmp(ti_fields{f_ix},'sample_rate')
            if trl_info.(ti_fields{f_ix})~=trl_infos{b_ix}.(ti_fields{f_ix}); error('sample_rate mismatch!'); end
        end
    end
end

%% Parameters
if ~isfield(proc_vars,'RT_std_thresh')
    proc_vars.RT_std_thresh = 3;
end
if ~isfield(proc_vars,'trial_lim_s')
    proc_vars.trial_lim_s = [-0.25 3];
end

%% Select channels and events of interest
if strcmp(proc_vars.event_type,'stim')
    events = trl_info.trl_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trl_info.rsp_onset;
else
    error(stract('ERROR: unknown event_type ',proc_vars.event_type));
end

% Convert trial_lim into samples
trial_lim = proc_vars.trial_lim_s*data.fsample;

%% Reject known artifacts
skip_rt  = find(trl_info.rt<0);     % no response (-1)
skip_bad = find(trl_info.hit<0);    % bad trial

% Convert visually bad epochs from full time to analysis_time
% Not necessary anymore since preproc epochs now include any kept from preclean
% %   NOTE: 1 keeps epochs that are only partially overlaping real data
% %   (i.e., will be trimmed to edges of analysis_time)
% bad_epochs = fn_convert_epochs_full2at(bad_epochs,SBJ_vars.analysis_time,...
%                                     strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean.mat'),1);
                                
% Toss epochs that overlap with bad_epochs from visual inspection
skip_vis = fn_find_trials_overlap_epochs(bad_epochs,1:size(data.trial{1},2),events,trial_lim);

% Find RT outliers
RT_mean = nanmean(trl_info.rt);
RT_std  = nanstd(trl_info.rt);
skip_rt_outlier = find(abs(trl_info.rt-RT_mean)>proc_vars.RT_std_thresh*RT_std);

% Check against RT bounds, toss late but only warn for early
RT_late = find(trl_info.rt>proc_vars.rt_bounds(2));
RT_early = find(trl_info.rt(trl_info.rt>0)<proc_vars.rt_bounds(1));
skip_rt_outlier = [skip_rt_outlier; RT_late; RT_early];

% Toss training trials
skip_training = find(trl_info.blk==0);

% Compile all a priori bad trials
skip_trial_ix = unique([skip_bad; skip_rt; skip_vis; skip_rt_outlier; skip_training]);
ok_trial_ix = setdiff(1:numel(trl_info.trl_n),skip_trial_ix);

%% Compile Bad Trials
% error('adjust trl_n for the training trials!');
trl_info_cln = trl_info;
trl_info_cln.event_type    = proc_vars.event_type;
trl_info_cln.trial_lim     = trial_lim;
trl_info_cln.trial_lim_s   = proc_vars.trial_lim_s;
trl_info_cln.RT_std_thresh = proc_vars.RT_std_thresh;

% Document bad trials
trl_info_cln.bad_trials.RT_bad = trl_info.trl_n(skip_rt);
trl_info_cln.bad_trials.RT_out = trl_info.trl_n(skip_rt_outlier);
trl_info_cln.bad_trials.visual = trl_info.trl_n(skip_vis);
trl_info_cln.bad_trials.bad    = trl_info.trl_n(skip_bad);
trl_info_cln.bad_trials.train  = trl_info.trl_n(skip_training);
trl_info_cln.bad_trials.all    = trl_info.trl_n(skip_trial_ix);

% Remove bad trials
n_trials = numel(trl_info.trl_n);
fields = fieldnames(trl_info_cln);
for f_ix = 1:numel(fields)
    if numel(trl_info_cln.(fields{f_ix}))==n_trials
        trl_info_cln.(fields{f_ix})(skip_trial_ix) = [];
    end
end

% Prevent formatting errors
trl_info_cln.rsp_onset = round(trl_info_cln.rsp_onset);

%% Print results
fprintf('==============================================================================================\n');
if ~isempty(RT_late)
    fprintf('WARNING! %i RTs > %f sec excluded!\n',numel(RT_late),proc_vars.rt_bounds(2));
end
if ~isempty(RT_early)
    fprintf('WARNING! %i RTs < %f sec excluded!\n',numel(RT_early),proc_vars.rt_bounds(1));
end
fprintf('Num trials excluded for training  : %i\n',length(skip_training));
fprintf('Num trials excluded for skip RT    : %i\n',length(skip_rt));
fprintf('Num trials excluded for outlier RT: %i\n',length(skip_rt_outlier));
fprintf('Num trials excluded by visual rej : %i\n',length(skip_vis));
fprintf('Num trials excluded for other     : %i\n',length(skip_bad));
fprintf('TOTAL TRIALS EXCLUDED A PRIORI    : %i\n',length(skip_trial_ix));
fprintf('TRIALS REMAINING: %i/%i\n',length(trl_info_cln.trl_n),length(trl_info.trl_n));
fprintf('==============================================================================================\n');

%% Save results
if save_it
    % Print results to logging file
    results_fname = [SBJ_vars.dirs.events SBJ '_behavior_rejection_results.txt'];
    r_file = fopen(results_fname,'a');
    fprintf(r_file,'%s\n',datestr(datetime));
    if ~isempty(RT_late)
        fprintf(r_file,'WARNING! %i RTs > %f sec excluded!\n',numel(RT_late),proc_vars.rt_bounds(2));
    end
    if ~isempty(RT_early)
        fprintf(r_file,'WARNING! %i RTs < %f sec excluded!\n',numel(RT_early),proc_vars.rt_bounds(1));
    end
    fprintf(r_file,'Num trials excluded for training  : %i\n',length(skip_training));
    fprintf(r_file,'Num trials excluded for skip RT   : %i\n',length(skip_rt));
    fprintf(r_file,'Num trials excluded for outlier RT: %i\n',length(skip_rt_outlier));
    fprintf(r_file,'Num trials excluded by visual rej : %i\n',length(skip_vis));
    fprintf(r_file,'Num trials excluded for other     : %i\n',length(skip_bad));
    fprintf(r_file,'TOTAL TRIALS EXCLUDED A PRIORI    : %i\n',length(skip_trial_ix));
    fprintf(r_file,'TRIALS REMAINING: %i/%i\n',length(trl_info_cln.trl_n),length(trl_info.trl_n));
    fclose(r_file);
    
    % Save the clean trl_info
    trl_info  = trl_info_cln;
    output_fname = [SBJ_vars.dirs.events SBJ '_trl_info_final.mat'];%',block_suffix,'.mat'];
    save(output_fname, '-v7.3', 'trl_info');
end

end
