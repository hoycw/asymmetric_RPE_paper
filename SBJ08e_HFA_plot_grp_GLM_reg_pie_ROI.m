function SBJ08e_HFA_plot_grp_GLM_reg_pie_ROI(SBJ_id, proc_id, stat_regs, hemi, roi_id,...
                                       plot_out, plt_id, save_fig, varargin)
%% Plot pie chart of overlap in mass GLM effects per ROI
%   COMMENTS NOT ADJSUTED SINCE STROOP COPY
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_regs [cell array] - {{'an_id1','model_id1','reg_lab1','stat_id1'},...
%                              {'an_idN','model_idN','reg_labN','stat_idN'}};
%       an_id [str] - analysis ID for preprocessing, filtering, etc.
%       model_id [str] - ID of the model to load
%       reg_lab [str] - label of the regressor (e.g., 'EV','sRPE','uRPE')
%       stat_id_n [str] - ID of statistical analysis
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   plt_id [str] - ID of the plotting variables (likely 'venn')
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?
%   save_fig [0/1] - save this figure?
%   fig_ftype [str] - file extension for figure saving

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Handle Variable Inputs & Defaults
% add %gm_thresh,z_thresh,plot_nsig,...
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        elseif strcmp(varargin{v},'atlas_id') && ischar(varargin{v+1})
            atlas_id = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');   fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('atlas_id','var');  atlas_id = 'Dx'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Define default options
eval(['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);

SBJs = fn_load_SBJ_list(SBJ_id);

% Organize IDs
if numel(stat_regs)==3; error('pie logic is not ready for venn with 3'); end
if numel(stat_regs) < 2 || numel(stat_regs) > 3; error('why venn?'); end
for st_ix = 1:numel(stat_regs)
    if numel(stat_regs{st_ix})~=4; error('need an_id, model_id, reg_lab, and stat_id in each stat_reg'); end
end
an_ids = cell(size(stat_regs)); stat_ids = cell(size(stat_regs));
model_ids = cell(size(stat_regs)); reg_labs = cell(size(stat_regs));
sts = cell(size(stat_regs));
for st_ix = 1:numel(stat_regs)
    an_ids{st_ix}    = stat_regs{st_ix}{1};
    model_ids{st_ix} = stat_regs{st_ix}{2};
    reg_labs{st_ix}   = stat_regs{st_ix}{3};
    stat_ids{st_ix}  = stat_regs{st_ix}{4};
    eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_ids{st_ix} '_vars.m']);
    sts{st_ix} = st;
    if ~strcmp(sts{st_ix}.an_style,'mGLM')
        error('This script is only for GLM analyses generating beta output structs!');
    end
    clear st
end
if all(strcmp(an_ids{1},an_ids))
    an_id_match = 1;
    an_id_cat = an_ids{1};
else
    an_id_match=0;
    an_id_cat = strjoin(an_ids,'-');
    warning('Not all an_ids match!');
end
if all(strcmp(model_ids{1},model_ids))
    model_id_match = 1;
    model_id_cat = model_ids{1};
else
    model_id_match=0;
    model_id_cat = strjoin(model_ids,'-');
    warning('Not all model_ids match!');
end
if all(strcmp(stat_ids{1},stat_ids))
    stat_id_match = 1;
    stat_id_cat = stat_ids{1};
else
    stat_id_match = 0;
    stat_id_cat = strjoin(stat_ids,'-');
    warning('Not all stat_ids match!');
end
reg_lab_cat = strjoin(reg_labs,'-');

% Venn colors
% if strcmp(reg_ids{1},'CNI') && strcmp(reg_ids{2},'pCNI') && strcmp(reg_ids{3},'PC')
%     color_id = 'ConCSPC';
% end
if exist('color_id','var')
    venn_colors = fn_venn_colors(numel(stat_regs), 'cond_id', color_id);
else
    venn_colors = fn_venn_colors(numel(stat_regs));
end
venn_colors = squeeze(venn_colors(logical(eye(size(venn_colors,1)))));
% venn_colors = cell(size(stat_conds));
% venn_colors{1} = plt.cmap(1,:);
% if numel(stat_conds)==3;
%     venn_colors{2} = plt.cmap(50,:);
%     venn_colors{3} = plt.cmap(end,:);
% else
%     venn_colors{2} = plt.cmap(end,:);
% end
[roi_list, roi_colors, roi_field] = fn_roi_label_styles(roi_id);

%% Prep report
% Create output dir
save_dir = [root_dir 'PRJ_Error/results/HFA/GRP/model_reg_pie/' stat_id_cat ...
            '/' model_id_cat '/' an_id_cat '/'];
if ~exist(save_dir,'dir'); [~] = mkdir(save_dir); end

% Create report
sig_report_fname = [save_dir SBJ_id '_' reg_lab_cat '_' atlas_id '_' roi_id '_' hemi '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
% Print IDs
title_str = repmat('%-40s',[1 numel(stat_regs)]);
fprintf(sig_report,[title_str '\n'],an_ids{:});
fprintf(sig_report,[title_str '\n'],model_ids{:});
fprintf(sig_report,[title_str '\n'],stat_ids{:});
fprintf(sig_report,[title_str '\n'],reg_labs{:});
fprintf(sig_report,[repmat('=',[1 40*numel(stat_regs)]) '\n']);
% Print result headings
fprintf(sig_report,[repmat('%-20s',[1 numel(stat_regs)+2]) '\n'],'label','ROI',reg_labs{:});
fprintf(sig_report,[repmat('-',[1 40*numel(stat_regs)]) '\n']);
result_str = ['%-20s%-20s' repmat('%-20i',[1 numel(stat_regs)]) '\n'];

%% Load Data
elec_sbj    = cell([numel(SBJs) 1]);
good_sbj    = true([numel(SBJs) 1]);
sig_roi_mat = cell([numel(SBJs) 1]);
roi_mat     = cell([numel(SBJs) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Prepare elec structs
    % Load elec struct
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    
    % Append SBJ name to labels
    orig_labels = elec_sbj{sbj_ix}.label;
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Check ROIs & hemisphere
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
    else
        % Check only hemisphere
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, [], []);
    end
    
    roi_mat{sbj_ix} = zeros([numel(elec_sbj{sbj_ix}.label) 1]);
    for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
        if any(strcmp(elec_sbj{sbj_ix}.label{ch_ix},roi_elecs))
            roi_mat{sbj_ix}(ch_ix) = find(strcmp(roi_list,elec_sbj{sbj_ix}.(roi_field){ch_ix}));
        end
    end
    
    %% Load Stats
    sig_mat = zeros([numel(elec_sbj{sbj_ix}.label) numel(stat_regs)]);
    if any(roi_mat{sbj_ix})
        for st_ix = 1:numel(stat_regs)
            % Load correct stat output structure, ensure match to elec
            if ~an_id_match || ~model_id_match || ~stat_id_match || st_ix==1
                if strcmp(reg_labs{st_ix},'actv')
                    error('not done actv yet');
                    % load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat'],'actv');
                elseif strcmp(reg_labs{st_ix},'RT')
                    error('RT not implemented');
                else
                    % Mass GLM
                    load([SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_ids{st_ix} '_' stat_ids{st_ix} '_' an_ids{st_ix} '.mat']);
                    if numel(elec_sbj{sbj_ix}.label)~=numel(beta.label) || ~all(strcmp(orig_labels,beta.label))
                        warning(['mismatched labels in elec and GLM beta struct for ' SBJ '!']);
                    end
                end
            end
            
            % Consolidate to binary sig/non-sig
            beta_reg_ix = strcmp(beta.feature,reg_labs{st_ix});
            for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
                beta_ch_ix = strcmp(beta.label,orig_labels{ch_ix});
                if any(squeeze(beta.qval(beta_reg_ix,beta_ch_ix,:))<=sts{st_ix}.alpha)
                    sig_mat(ch_ix,st_ix) = 1;
                end
            end
            if ~an_id_match || ~model_id_match || ~stat_id_match || st_ix==numel(stat_regs); clear beta; end
            clear beta_reg_ix
        end
        
        %% Keep intersection of significant and ROI matched electrodes
        if any(any(sig_mat(roi_mat{sbj_ix}~=0,:),2))
            sig_idx = any(sig_mat,2);
            roi_idx = roi_mat{sbj_ix}~=0;
            sig_roi_mat{sbj_ix} = sig_mat(roi_idx & sig_idx,:).*roi_mat{sbj_ix}(roi_idx & sig_idx);
            fprintf('\t%s has %i sig channels in %s hemi %s\n',SBJ,size(sig_roi_mat{sbj_ix},1),atlas_id,hemi);
        else
            good_sbj(sbj_ix) = false;
            fprintf(2,'\t%s has %i channels in %s hemi %s, but none are significant\n',...
                SBJ,sum(roi_mat{sbj_ix}~=0),atlas_id,hemi);
        end
    else
        good_sbj(sbj_ix) = false;
        fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJ,atlas_id,hemi);
    end
    %% Report all elecs, even non-ROI and non-sig
    for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
        fprintf(sig_report,result_str,elec_sbj{sbj_ix}.label{ch_ix},...
            elec_sbj{sbj_ix}.(roi_field){ch_ix},sig_mat(ch_ix,:));
    end
    clear SBJ SBJ_vars SBJ_vars_cmd sig_elecs sig_ix cond_ix sig_mat roi_idx sig_idx
end

% Finish report and save sig_mat
fclose(sig_report);
% sig_elecs_fname = [save_dir SBJ_id '_' atlas_id '_' roi_id '_sig_elecs.mat'];
% save(sig_elecs_fname,'-v7.3','stat_conds','SBJs','elec_sbj','sig_mat');

%% Combine elec structs
% elec = ft_appendsens([],elec_sbj{good_sbj});
% elec.roi_color = vertcat(all_roi_colors{:});    % appendsens strips that field
% elec.(roi_field)       = all_roi_labels{cond_ix};    % appendsens strips that field

sig_roi_all = vertcat(sig_roi_mat{:});

%% Compute Stat Overlaps
% Find number of comparisons
pairs = nchoosek(1:numel(stat_regs),2);
n_groups = numel(stat_ids)+size(pairs,1);
if numel(stat_regs)==3; n_groups = n_groups+1; end

sig_cnt     = zeros([n_groups 1]);                  % # significant elecs per group
sig_type    = zeros([n_groups 1]);                  % individual (1), pair (2), or triple (3) effect
sig_grps    = cell([n_groups 1]);                   % effect assignments per group
sig_cnt_roi = zeros([n_groups numel(roi_list)]);    % count per ROI per group
pie_legend  = cell(size(stat_regs));
roi_legends = cell([numel(stat_regs) numel(roi_list)]);
% Add main effects
grp_ix = 0;
for main_ix = 1:numel(stat_regs)
    grp_ix = grp_ix + 1;
    sig_type(grp_ix) = 1;
    non_main_ix = setdiff(1:numel(stat_ids),main_ix);
    sig_cnt(grp_ix) = sum(sig_roi_all(:,main_ix)~=0 & sig_roi_all(:,non_main_ix)==0);
    sig_grps{grp_ix} = main_ix;
    pie_legend{grp_ix} = [reg_labs{main_ix} ' = ' num2str(sig_cnt(grp_ix))];
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(grp_ix,roi_ix) = sum(sig_roi_all(:,main_ix)==roi_ix & sig_roi_all(:,non_main_ix)~=roi_ix);
        roi_legends{grp_ix,roi_ix} = [reg_labs{main_ix} ' = ' num2str(sig_cnt_roi(grp_ix,roi_ix))];
    end
end

% Add pairs
for p_ix = 1:size(pairs,1)
    grp_ix = grp_ix + 1;
    sig_type(grp_ix) = 2;
    sig_grps{grp_ix} = pairs(p_ix,:);
    sig_cnt(grp_ix) = sum(sig_roi_all(:,pairs(p_ix,1))~=0 & sig_roi_all(:,pairs(p_ix,2))~=0);
    pie_legend{grp_ix} = [reg_labs{pairs(p_ix,1)} '(#' num2str(pairs(p_ix,1)) ') + '...
        reg_labs{pairs(p_ix,2)} '(#' num2str(pairs(p_ix,2)) ') = ' num2str(sig_cnt(grp_ix))];
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(grp_ix,roi_ix) = sum(sig_roi_all(:,pairs(p_ix,1))==roi_ix & ...
                                  sig_roi_all(:,pairs(p_ix,2))==roi_ix);
        roi_legends{grp_ix,roi_ix} = [reg_labs{pairs(p_ix,1)} '(#' num2str(pairs(p_ix,1)) ') + '...
                            reg_labs{pairs(p_ix,2)} '(#' num2str(pairs(p_ix,2)) ') = ' num2str(sig_cnt_roi(grp_ix,roi_ix))];
    end
end

% Add triples
% if numel(stat_regs)==3
%     sig_cnt(end) = sum(all(sig_roi_all,2));
%     sig_type(end) = 3;
%     sig_grps{end} = 1:3;
% %     venn_legend{end} = ['ALL (n=' num2str(sig_cnt(end)) ')'];
%     for roi_ix = 1:numel(roi_list)
%         sig_cnt_roi(end,roi_ix) = sum(all(sig_roi_all==roi_ix,2));
%     end
% end

%% Plot Overall Pie Chart
fig_name = [SBJ_id '_pie_' atlas_id '_' roi_id '_' stat_id_cat];
if numel(stat_regs)==3; pos = plt.fig_pos{2};
else pos = plt.fig_pos{1};
end
f = figure('Name',fig_name,'units','normalized','outerposition',pos);

pie_obj = pie(sig_cnt,pie_legend);
for pie_ix = 1:numel(sig_cnt)
    pie_obj(pie_ix*2).FontSize = plt.text_sz;
end
% grp_ix = 0;
% for main_ix = 1:numel(stat_regs)
%     grp_ix = grp_ix + 1;
% %     pie_obj(grp_ix*2-1).BackgroundColor = venn_colors{grp_ix};
%     pie_obj(grp_ix*2).FontSize = plt.text_sz;
% end
% for p_ix = 1:size(pairs,1)
%     grp_ix = grp_ix + 1;
%     pie_obj(grp_ix*2-1).BackgroundColor = venn_colors{grp_ix};
%     pie_obj(grp_ix*2).FontSize = plt.text_sz;
% end
title([atlas_id '-' roi_id ': # SBJ = ' num2str(numel(SBJs)) ', # Sig+ROI Elecs = '...
       num2str(size(sig_roi_all,1)) ' / ' num2str(size(vertcat(roi_mat{:}),1)) ...
       ' (' num2str(size(sig_roi_all,1)/size(vertcat(roi_mat{:}),1),'%.4f') ')'],...
       'FontSize',plt.title_sz);

if save_fig
    fig_fname = [save_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Venn Diagrams per ROI
for roi_ix = 1:numel(roi_list)
    fig_name = [SBJ_id '_pie_' atlas_id '_' roi_id '_' stat_id_cat '_' roi_list{roi_ix}];
    f = figure('Name',fig_name,'units','normalized','outerposition',pos);
    
    pie_obj = pie(sig_cnt_roi(:,roi_ix),roi_legends(:,roi_ix));
    for p_ix = 1:numel(sig_cnt)
        %     pie(p_ix*2-1).BackgroundColor = venn_colors{1};
        pie_obj(p_ix*2).FontSize = plt.text_sz;
    end
    
    title([atlas_id '-' roi_id '-' roi_list{roi_ix} ': # SBJ = ' num2str(numel(SBJs)) ...
        ', # Sig Elecs = ' num2str(sum(any(sig_roi_all==roi_ix,2))) ' / ' num2str(sum(vertcat(roi_mat{:})==roi_ix)) ...
       ' (' num2str(sum(any(sig_roi_all==roi_ix,2))/sum(vertcat(roi_mat{:})==roi_ix),'%.4f') ')'],...
       'FontSize',plt.title_sz,'Color',roi_colors{roi_ix});
    
    if save_fig
        fig_fname = [save_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
