function SBJ08l_reversal_points(SBJs, an_id, model_id, stat_id, chpval_type)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname)
wvals = [];
for sj = 1:length(SBJs)
    SBJ = SBJs{sj};
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
    load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
    load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'.mat']);
    [roi_list,~, roi_field] = fn_roi_label_styles(roi_id);
    load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    [~,EV_order] = sort(model(:,1));
    r_ix = [];
    for r = 1:numel(roi_list)
        cr_ix = find(strcmp(elec.(roi_field), roi_list{r}));
        r_ix = [r_ix; cr_ix];
    end
    hfa.label = hfa.label(r_ix);
    hfa.powspctrm = hfa.powspctrm(EV_order,r_ix,:,:);
    for ch = 1:numel(hfa.label)
        cchan = [SBJ ' ' hfa.label{ch}];
        chix = [0,0];
        for r = 1:numel(beta_chan.label)
            if ismember(cchan, beta_chan.chan_label{r})
                chix(1) = r;
                chix(2) = find(strcmp(beta_chan.chan_label{r}, cchan));
            end
        end
        sig = sum(beta_chan.pvals{chix(1)}(chix(2),3:4,:) < .05,'all') > 0;
        tidx = hfa.time >= 0.2 & hfa.time <= 0.45; 
        if sig > 0
           wvals(end+1) = sum(model(EV_order,1) .* abs(squeeze(mean(hfa.powspctrm(:,ch,:,tidx),4)))) / size(model,1);
          
            figure;
            plot(model(EV_order,1),squeeze(mean(hfa.powspctrm(:,ch,:,tidx),4)),'.')
            title(cchan)         
        end
    end
    clear dirs_fields field_ix SBJ_id SBJ_vars hfa hfa_win
end
