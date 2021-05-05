%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%% Load SBJ list
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);

%% ERP Analysis and Plotting
an_id   = 'ERP_F_trl2to1_stat1';
fig_vis = 'off';
for s = 1:numel(SBJs)
%     conditions = 'EzOut';
%     SBJ06a_ERP_stats(SBJs{s},conditions,proc_id,an_id);
%     conditions = 'HdOut';
%     SBJ06a_ERP_stats(SBJs{s},conditions,proc_id,an_id);
%     
%     conditions = 'DifOut';
%     plt_id     = 'ts_trl2to1000_errbr';
%     SBJ06b_ERP_plot_stats(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,fig_vis)
%     close all;
end

%% Compute High Frequency Activity via Sun Grid Engine
% HFA filtering (especially with multitapers) take a long time, so use SGE
% On the cluster, run:
%   `cd /${root_dir}/PRJ_Error/scripts/
%   `submit -s SGE/calls/SBJ07a_HFA_save.sh -f SBJ_lists/list.sbj \
%       -o SGE/opts/SBJ07a_save_HGm_F25t121_zbtS_sm0_l1_wn100.opts`
% INPUTS:
%   -s script.sh: shell script to cd to scripts, load the SBJ_list and options, 
%       create a one line .m script to call the function, open MATLAB, run
%       the script, and delete it before exiting the job
%   -f SBJ_list.sbj: text file with SBJ names, one per line, will run one
%       qsub job per SBJ
%   -o file.opts: options for SGE job
%       -o /full/path/to/standard_out output file
%       -e /full/path/to/standard_error output file
%       -v bash script ${variables} that can be used within script.sh (e.g., an_id, proc_id, etc.)
%       -l mem_free=16G - request memory for a node
%       -N name of the job (e.g., to view progress via `qstat`)

%% Test and Plot Active HFA Channels
conditions = 'DifOut';
proc_id    = 'main_ft';
an_id      = 'HGm_F25t121_zbtS_sm0_l1_wn100';
actv_win   = 100;
atlas_id   = 'Dx';

plt_id     = 'stack_F2to12_evnt_c5';%'stack_S2to3_evnt_c5';
save_fig   = 1;
fig_ftype  = 'png';
fig_vis    = 'on';

for s = 1:numel(SBJs)
    % Compute and Save Active Channels with High Frequency Activity
%     SBJ07b_HFA_actv(SBJs{s},proc_id,an_id,actv_win);
    
    % Plot Single Trial HFA Stacks
    SBJ07c_HFA_plot_stack_mean(SBJs{s}, conditions, proc_id, an_id, actv_win, plt_id, save_fig,...
        'atlas_id', atlas_id, 'fig_vis', fig_vis, 'fig_ftype', fig_ftype);
    close all;
end

