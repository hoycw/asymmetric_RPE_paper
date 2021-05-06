function SBJ07b_HFA_actv(SBJ,proc_id,an_id,actv_win)
%% Calculates activation relative to baseline via point wise t-test with FDR
% INPUTS:
%   SBJ [str] - ID of subject
%   proc_id [str] - name of processing pipeline (proc_vars)
%   an_id [str] - ID of the analysis parameters to use
%   actv_win [int] - length in ms of minimum difference from baseline to be active
% OUTPUTS:
%   actv [struct]:
%       .an_id [str] - ID of analysis used
%       .actv_win [int] - length in ms of minimum difference from baseline to be active
%       .pvals [float] - [n_chan, n_time] matrix of p values
%       .qvals [float] - [n_chan, n_time] matrix of p values corrected via FDR
%       .actv_ch [0/1] - binary indicator is any channel is active
%       .actv_epochs [cell array] - [n_epochs, 2] array per channel with 
%           [start, stop] of active epochs in seconds

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
if ischar(actv_win); actv_win = str2num(actv_win); end

%% Set up paths
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
% eval(['run /home/knight/hoycw/PRJ_Error/scripts/proc_vars/' proc_id '_vars.m']);

% Load Data
hfa_fname = [SBJ_vars.dirs.proc SBJ '_ROI_' proc_id '_' an_id '.mat'];
load(hfa_fname);
load([SBJ_vars.dirs.events SBJ '_bhv_' proc_id '_final.mat']);

% Convert
sample_rate = (numel(hfa.time)-1)/(hfa.time(end)-hfa.time(1));
actv_win_samples = sample_rate*actv_win/1000;

%% Run Statistics
fprintf('==================================================================\n');
fprintf('--------------------- Computing Active Channels ------------------\n');
fprintf('==================================================================\n');

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = an.trial_lim_s;
hfa_stat = ft_selectdata(cfg_trim,hfa);

actv = struct();
actv.an_id       = an_id;
actv.actv_win    = actv_win;
actv.pvals       = nan([numel(hfa_stat.label) numel(hfa_stat.time)]);
actv.qvals       = nan([numel(hfa_stat.label) numel(hfa_stat.time)]);
actv.actv_ch     = zeros(size(hfa_stat.label));
actv.actv_epochs = cell(size(hfa_stat.label));
for ch_ix = 1:numel(hfa_stat.label)
    % Compute t-test per time point
    for time_ix = 1:numel(hfa_stat.time)
        [~, actv.pvals(ch_ix,time_ix)] = ttest(squeeze(hfa_stat.powspctrm(:,ch_ix,1,time_ix)));
    end
    
    % Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch (e.g., ROF8-9
%     in IR32), so I'm trying the below function
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [~, ~, ~, actv.qvals(ch_ix,:)] = fdr_bh(actv.pvals(ch_ix,:));
    actv_mask = actv.qvals(ch_ix,:)<=0.05;
    actv_chunks = fn_find_chunks(actv_mask);
    actv_chunks(actv_mask(actv_chunks(:,1))==0,:) = [];
    actv_chunk_len = diff(actv_chunks,1,2)+1;
    actv_epochs = actv_chunks(actv_chunk_len > actv_win_samples,:);
    if ~isempty(actv_epochs)
        actv.actv_ch(ch_ix) = 1;
        actv.actv_epochs{ch_ix} = hfa_stat.time(actv_epochs);
    end
end

%% Save Results
out_fname = strcat(hfa_fname(1:end-4),'_actv_mn',num2str(actv_win),'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',out_fname);
fprintf('===================================================\n');
save(out_fname,'-v7.3','actv');

end
