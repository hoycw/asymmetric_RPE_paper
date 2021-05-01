function SBJ07ab_HFA_actv(SBJ,an_id,actv_win)
% Calculates activation relative to baseline via point wise t-test with FDR

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
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
load(hfa_fname);
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));

%% Run Statistics
fprintf('===================================================\n');
fprintf('--------------------- Statistics ------------------\n');
fprintf('===================================================\n');

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = stat_lim;
hfa_stat = ft_selectdata(cfg_trim,hfa);

actv_ch = {};
actv_ch_epochs = {};
for ch_ix = 1:numel(hfa_stat.label)
    % Compute t-test per time point
    n_tests   = size(hfa_stat.powspctrm,4);
    pvals     = NaN([1 n_tests]);
    for time_ix = 1:n_tests
        [~, pvals(time_ix)] = ttest(squeeze(hfa_stat.powspctrm(:,ch_ix,1,time_ix)));
    end
    
    % Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch (e.g., ROF8-9
%     in IR32), so I'm trying the below function
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [~, ~, ~, qvals] = fdr_bh(pvals);
    actv_mask = qvals<=0.05;
    actv_chunks = fn_find_chunks(actv_mask);
    actv_chunks(actv_mask(actv_chunks(:,1))==0,:) = [];
    actv_chunk_sz = diff(actv_chunks,1,2)+1;
    actv_epochs = actv_chunks(actv_chunk_sz > actv_win,:);
    if ~isempty(actv_epochs)
        actv_ch = {actv_ch{:} hfa_stat.label{ch_ix}};
        actv_ch_epochs = {actv_ch_epochs{:}, hfa_stat.time(actv_epochs)};
    end
end

%% Save Results
out_fname = strcat(hfa_fname(1:end-4),'_actv_mn',num2str(actv_win),'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',out_fname);
fprintf('===================================================\n');
save(out_fname,'-v7.3','actv_ch','actv_ch_epochs','pvals','qvals');

end
