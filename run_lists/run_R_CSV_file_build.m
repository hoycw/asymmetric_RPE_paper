%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Full SBJ list (completed and ready to run)
SBJs      = {'CP24','IR57','IR68'};
proc_id   = 'main_ft';
an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_ids  = {'DifOutDO_F02t05_WL03_WS01','DifOutDO_F02t05_WL01_WS01','DifOutDO_F01t06_WL01_WS01','DifOutDO_F01t06_WL005_WS005'};
atlas_id  = 'Dx';
roi_id    = 'main3';

%% Prepare for inspections
for s = 1:numel(SBJs)
    fprintf('============= SBJ: %s ===============\n',SBJs{s});
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    for st_ix = 1:numel(stat_ids)
        fprintf('============= Stat: %s ===============\n',stat_ids{st_ix});
        % Build design predictors
        [design,des_names,des_vals] = fn_build_design_predictors(SBJs{s},stat_ids{st_ix});
        
        %     % Build ROI predictors
        %     [roi_info,roi_names,roi_vals] = fn_build_roi_predictors(SBJs{s},atlas_id,roi_id);
        
        % Build data feature matrix
        [data,data_names,data_vals] = fn_build_data_features(SBJs{s},proc_id,an_id,stat_ids{st_ix},atlas_id,roi_id);
        
        % Combine matrices
        full = horzcat(design,data);
        full_names = horzcat(des_names,data_names);
        full_vals  = horzcat(des_vals,data_vals);
        matrix = vertcat(full_names,num2cell(full));
        
        % Convert NaN to 'NA' for R processing
        nan_idx = isnan(full);
        for r_ix = 1:size(full,1)
            if any(nan_idx(r_ix,:))
                nan_ix = find(nan_idx(r_ix,:));
                for ix = find(nan_idx(r_ix,:))
                    matrix{r_ix+1,ix} = 'NA';
                end
            end
        end
        
        % Save CSV file
        csv_fname = [SBJ_vars.dirs.proc SBJs{s} '_ROI_' stat_ids{st_ix} '_' an_id '_' atlas_id '_' roi_id '_full_model.csv'];
        fprintf('============= Saving %s ===============\n',csv_fname);
        % Write header
        fid = fopen(csv_fname,'w');
        header = strjoin(full_names,',');
        fprintf(fid,'%s\n',header);
        fclose(fid);
        
        % Append data to end of file
        dlmwrite(csv_fname,full,'-append');
        fprintf('\n\n\n');
    end
end



