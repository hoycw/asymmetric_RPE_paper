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
conditions = 'EzOut';%'DifFB';
proc_id    = 'main_ft';
an_id      = 'ERP_F25t1';
atlas_id   = 'Dx';

plt_id     = 'stack_F2t1_evnt_c5';
save_fig   = 1;
fig_ftype  = 'png';
fig_vis    = 'off';

for s = 1:numel(SBJs)
    SBJ06a_ERP_save(SBJs{s},proc_id,an_id);
    
    SBJ06b_ERP_plot_stack_mean(SBJs{s}, conditions, proc_id, an_id, plt_id, save_fig,...
        'atlas_id', atlas_id, 'fig_vis', fig_vis, 'fig_ftype', fig_ftype);
    close all;
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
conditions = 'DifFB';
proc_id    = 'main_ft';
% an_id      = 'HGh_F25t121_zbtS_sm0_l0';%'HGh_S25t301_zbtS_sm0_l1';%'HGh_F25t121_zbtS_sm0_l1';%
an_id      = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGm_S25t301_zbtS_sm0_l1_wn100';%
actv_win   = 100;
atlas_id   = 'Dx';

plt_id     = 'stack_F2t1_evnt_c5';%'stack_S2t3_evnt_c5';%
save_fig   = 1;
fig_ftype  = 'png';
fig_vis    = 'on';

for s = 9%1:numel(SBJs)
    % Compute and Save Active Channels with High Frequency Activity
%     SBJ07b_HFA_actv(SBJs{s},proc_id,an_id,actv_win);
    
    % Plot Single Trial HFA Stacks
    SBJ07c_HFA_plot_stack_mean(SBJs{s}, conditions, proc_id, an_id, actv_win, plt_id, save_fig,...
        'atlas_id', atlas_id, 'fig_vis', fig_vis, 'fig_ftype', fig_ftype, 'elec_lab', {'LIN3-4'});
    %close all;
end

