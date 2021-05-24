%% Print number of electrodes per ROI and per SBJ
SBJ_id = 'exit';
SBJs = fn_load_SBJ_list(SBJ_id);

roi_id = 'MPFCINS';
atlas_id = 'Dx';
proc_id = 'main_ft';

[roi_list, ~, roi_field] = fn_roi_label_styles(roi_id);
roi_cnt  = zeros([numel(roi_list) numel(SBJs)]);
for s = 1:numel(SBJs)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    elec_fname = [SBJ_vars.dirs.recon,SBJs{s},'_elec_',proc_id,'_pat_',atlas_id,'_final.mat'];
    load(elec_fname);
    for roi_ix = 1:numel(roi_list)
        roi_cnt(roi_ix,s) = sum(strcmp(elec.(roi_field),roi_list{roi_ix}));
    end
end

for roi_ix = 1:numel(roi_list)
    fprintf('%s total = %d; per SBJ = %.2f (%d - %d); nSBJ = %d\n',roi_list{roi_ix},...
        sum(roi_cnt(roi_ix,:)),sum(roi_cnt(roi_ix,:))/numel(SBJs),...
        min(roi_cnt(roi_ix,:)),max(roi_cnt(roi_ix,:)),sum(roi_cnt(roi_ix,:)~=0));
end

fprintf('elec/SBJ = %.2f +/- %.2f\n',mean(sum(roi_cnt,1)),std(sum(roi_cnt,1)));