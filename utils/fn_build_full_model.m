function [model,model_names] = fn_build_full_model(SBJ,design,design_names,proc_id,an_id,stat_id,atlas_id,roi_id)
%% Build design matrix of predictors for regression analysis
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Load Data
fprintf('================== Loading Data =======================\n');
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat'));

% Check if more than one frequency, error for now
if numel(hfa.freq)>1
    error('HFA has more than one frequency, can''t run on that for now!');
end

%% Load Elec and ROI info
fprintf('================== ROI Overlap =======================\n');
load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);

% Load all ROI info
[roi_list, ~] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

% Find ROI matched elecs
roi_idx = false([numel(elec.label) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
    roi_idx(:,roi_ix) = strcmp(elec.(roi_field),roi_list{roi_ix});
end
roi_idx = any(roi_idx,2);

%% Select data in stat window
cfg_trim = [];
cfg_trim.channel = intersect(elec.label(roi_idx),hfa.label);
elec = fn_select_elec(cfg_trim,elec);
cfg_trim.latency = [st.stat_lim(1) st.stat_lim(2)+0.001];
hfa = ft_selectdata(cfg_trim,hfa);

if ~all(strcmp(elec.label,hfa.label)), error('mismatch hfa elec label order'); end

%% Log transform HFA


%% Average HFA in Sliding Windows
fprintf('================== Averaging HFA within Windows =======================\n');
% Sliding window parameters
win_lim    = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),...
    round(st.win_len*trl_info.sample_rate),...
    round(st.win_step*trl_info.sample_rate));
win_center = round(mean(win_lim,2));

% Average in windows
hfa_win = zeros([size(hfa.powspctrm,1) size(hfa.powspctrm,2) size(win_lim,1)]);%size(hfa.powspctrm,3)
win_names = cell(size(win_center));
for w_ix = 1:size(win_lim,1)
    hfa_win(:,:,w_ix) = squeeze(nanmean(hfa.powspctrm(:,:,1,win_lim(w_ix,1):win_lim(w_ix,2)),4));
    win_names{w_ix} = num2str(hfa.time(win_center(w_ix)),'%.2f');
end

%% Raster out into column format
model_names = [design_names {'channel', 'ROI', 'Time', 'HFA'}];
trl_n_ix  = find(strcmp(design_names,'trl_n'));
model     = cell(size(model_names));
row_ix   = 0;
for roi_ix = 1:numel(roi_list)
    ch_idx = find(strcmp(elec.(roi_field),roi_list{roi_ix}));
    for ch_ix = 1:numel(ch_idx)
        for time_ix = 1:size(win_lim,1)
            for trl_ix = 1:numel(trl_info.trl_n)
                row_ix = row_ix + 1;
                design_ix = find([design{:,trl_n_ix}]==trl_info.trl_n(trl_ix));
                model(row_ix,:) = [design(design_ix,:) elec.label(ch_idx(ch_ix)) roi_list(roi_ix) win_names(time_ix) hfa_win(trl_ix,ch_idx(ch_ix),time_ix)];
            end
        end
    end
end        

end