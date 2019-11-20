%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Full SBJ list (completed and ready to run)
SBJs      = {'CP24','IR57','IR68'};%'CP24'};%{
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_ids  = {'DifOutDO_F02t05_WL03_WS01','DifOutDO_F01t06_WL005_WS005'};%'DifOutDO_F02t05_WL01_WS01','DifOutDO_F01t06_WL01_WS01',
atlas_id  = 'Dx';
roi_id    = 'main3';

log_transform = 1;

%% Prepare for inspections
for s = 1:numel(SBJs)
    fprintf('============= SBJ: %s ===============\n',SBJs{s});
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    for st_ix = 1:numel(stat_ids)
        fprintf('============= Stat: %s ===============\n',stat_ids{st_ix});
        % Build design predictors
        [design,design_names] = fn_build_design_predictors(SBJs{s},stat_ids{st_ix});
        % Help steve with naming
        for name_ix = 1:numel(design_names)
            design_names{name_ix} = strrep(design_names{name_ix},'Dif*Out','Surp');
        end

        % Build data feature matrix
        [model,model_names] = fn_build_full_model(SBJs{s},design,design_names,...
                                            proc_id,an_id,stat_ids{st_ix},atlas_id,roi_id,log_transform);
        
        % Save CSV file
        csv_fname = [SBJ_vars.dirs.proc SBJs{s} '_ROI_' stat_ids{st_ix} '_' an_id '_' atlas_id '_' roi_id '_full_model.csv'];
        fprintf('============= Saving %s ===============\n',csv_fname);
%         % Write header
%         fid = fopen(csv_fname,'w');
%         header = strjoin(model_names,',');
%         fprintf(fid,'%s\n',header);
%         fclose(fid);
%         dlmwrite(csv_fname,model,'-append');
        
        % Append data to end of file
        % Convert cell to a table and use first row as variable names
        table = cell2table(model,'VariableNames',model_names);
        % Write the table to a CSV file
        writetable(table,csv_fname)
        fprintf('\n\n\n');
    end
end



