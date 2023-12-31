function SBJ08e_HFA_plot_grp_GLM_reg_cnfmat_ROI(SBJ_id, proc_id, stat_regs, hemi, roi_id,...
                                       plot_out, plt_id, save_fig, varargin)
%% Plot confusion matrix of overlap in mass GLM effects per ROI
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
addpath([app_dir 'cbrewer/']);
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
if numel(stat_regs) < 2 || numel(stat_regs) > 3; error('why venn?'); end
for st_ix = 1:numel(stat_regs)
    if numel(stat_regs{st_ix})~=4; error('need an_id, model_id, reg_lab, and stat_id in each stat_reg'); end
end
an_ids = cell(size(stat_regs)); stat_ids = cell(size(stat_regs));
model_ids = cell(size(stat_regs)); reg_labs = cell(size(stat_regs));
sts = cell(size(stat_regs)); reg_ixs = nan(size(stat_regs));
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
save_dir = [root_dir 'PRJ_Error/results/HFA/model_reg_cnfmat/' stat_id_cat ...
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
roi_all     = vertcat(roi_mat{:});

%% Compute Stat Overlaps
% Find number of comparisons
pairs = nchoosek(1:numel(stat_regs),2);

% Add main effects
sig_cnt          = nan(numel(stat_regs));                  % # significant elecs per group
sig_cnt_norm     = nan(numel(stat_regs));                  % # significant elecs per group
sig_cnt_roi      = nan([numel(stat_regs) numel(stat_regs) numel(roi_list)]);    % count per ROI per group
sig_cnt_norm_roi = nan([numel(stat_regs) numel(stat_regs) numel(roi_list)]);    % count per ROI per group
for main_ix = 1:numel(stat_regs)
    sig_cnt(main_ix,main_ix) = sum(sig_roi_all(:,main_ix)~=0);
    sig_cnt_norm(main_ix,main_ix) = sum(sig_roi_all(:,main_ix)~=0)/sum(roi_all~=0);
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(main_ix,main_ix,roi_ix) = sum(sig_roi_all(:,main_ix)==roi_ix);
        sig_cnt_norm_roi(main_ix,main_ix,roi_ix) = sum(sig_roi_all(:,main_ix)==roi_ix)/sum(roi_all==roi_ix);
    end
end

% Add pairs
for p_ix = 1:size(pairs,1)
    sig_cnt(pairs(p_ix,1),pairs(p_ix,2)) = sum(sig_roi_all(:,pairs(p_ix,1))~=0 & sig_roi_all(:,pairs(p_ix,2))~=0);
    sig_cnt_norm(pairs(p_ix,1),pairs(p_ix,2)) = ...
        sum(sig_roi_all(:,pairs(p_ix,1))~=0 & sig_roi_all(:,pairs(p_ix,2))~=0)/sum(roi_all~=0);
    for roi_ix = 1:numel(roi_list)
        sig_cnt_roi(pairs(p_ix,1),pairs(p_ix,2),roi_ix) = sum(sig_roi_all(:,pairs(p_ix,1))==roi_ix & ...
                                  sig_roi_all(:,pairs(p_ix,2))==roi_ix);
        sig_cnt_norm_roi(pairs(p_ix,1),pairs(p_ix,2),roi_ix) = sum(sig_roi_all(:,pairs(p_ix,1))==roi_ix & ...
                                  sig_roi_all(:,pairs(p_ix,2))==roi_ix)/sum(roi_all==roi_ix);
    end
end

%% Plot Overall Confusion Matrix Diagram
fig_name = [SBJ_id '_cnfmat_' atlas_id '_' roi_id '_' stat_id_cat];
f = figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.5 0.5]);

% Plot Full Count
subplot(1,2,1);
cmap = cbrewer('seq','Reds',max(sig_cnt(:))-min(sig_cnt(:)));
imagesc(sig_cnt,'AlphaData',~isnan(sig_cnt));
cbar = colorbar;
cbar.Limits = [0 max(sig_cnt(:))];
colormap(cmap);
set(gca,'YDir','normal');
for reg_ix1 = 1:numel(stat_regs)
    for reg_ix2 = reg_ix1:numel(stat_regs)
            text(reg_ix2,reg_ix1,num2str(sig_cnt(reg_ix1,reg_ix2)),...
                'FontSize',16,'HorizontalAlignment','center');
    end
end
xticks(1:numel(stat_regs));
xticklabels(reg_labs);
yticks(1:numel(stat_regs));
yticklabels(reg_labs);
title(['Sig. Electrode Counts in ' roi_id ' ROIs']);
set(gca,'FontSize',16);

% Plot Normalized Count
subplot(1,2,2);
cmap = cbrewer('seq','Reds',max(sig_cnt_norm(:))-min(sig_cnt_norm(:)));
imagesc(sig_cnt_norm,'AlphaData',~isnan(sig_cnt_norm));
cbar = colorbar;
cbar.Limits = [0 max(sig_cnt_norm(:))];
colormap(cmap);
set(gca,'YDir','normal');
for reg_ix1 = 1:numel(stat_regs)
    for reg_ix2 = reg_ix1:numel(stat_regs)
            text(reg_ix2,reg_ix1,num2str(sig_cnt(reg_ix1,reg_ix2)),...
                'FontSize',16,'HorizontalAlignment','center');
    end
end
xticks(1:numel(stat_regs));
xticklabels(reg_labs);
yticks(1:numel(stat_regs));
yticklabels(reg_labs);
title(['% Sig. Electrodes in ' roi_id ' ROIs']);
set(gca,'FontSize',16);

if save_fig
    fig_fname = [save_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Venn Diagrams per ROI
for roi_ix = 1:numel(roi_list)
    fig_name = [SBJ_id '_venn_' atlas_id '_' roi_id '_' stat_id_cat '_' roi_list{roi_ix}];
    f = figure('Name',fig_name,'units','normalized','outerposition',pos);
    ax = gca; hold on;
    
    if ~any(sig_cnt_roi(:,roi_ix)==0)
        % f = venn(Z), where:
        %   n=2: [z1 z2 z12]
        %   n=3: [z1 z2 z3 z12 z13 z23 z123]
        [v, v_stats] = venn(sig_cnt_roi(:,roi_ix),'FaceColor',venn_colors);
        text_spacers = zeros(size(v_stats.ZoneCentroid));
    else
        % Not all circles overlap, venn fails
        % Plot circles for each method
        [~,max_idx] = sort(sig_cnt_roi(:,roi_ix),'descend');
        scat_space = sum(sig_cnt_roi(max_idx(1:2),roi_ix))*cos(pi/4);
        scat_locs  = [-scat_space 0; scat_space 0; 0 scat_space];
        v = gobjects(size(stat_regs));
        v_stats.ZoneCentroid = zeros([n_groups 2]);
        for st_ix = 1:numel(stat_regs)
            v(st_ix) = patch(sig_cnt_roi(st_ix,roi_ix)*cos(plt.ang)+scat_locs(st_ix,1),...
                             sig_cnt_roi(st_ix,roi_ix)*sin(plt.ang)+scat_locs(st_ix,2),...
                             venn_colors{st_ix});
            set(v(st_ix),'FaceAlpha',plt.alpha);
            v_stats.ZoneCentroid(st_ix,:) = scat_locs(st_ix,:);
        end
        ax.XLim = [-scat_space-max(sig_cnt_roi(:,roi_ix))-plt.ax_fudge ...
                    scat_space+max(sig_cnt_roi(:,roi_ix))+plt.ax_fudge];
        ax.YLim = [-scat_space-max(sig_cnt_roi(:,roi_ix))-plt.ax_fudge ...
                    scat_space+max(sig_cnt_roi(:,roi_ix))+plt.ax_fudge];
        
        % Plot lines to show overlap
        grp_idx = find(sig_type==2);
        for p_ix = 1:size(pairs,1)
            grp_ix = grp_idx(p_ix);
            if sig_cnt_roi(grp_ix,roi_ix)~=0
                xs = [scat_locs(sig_grps{grp_ix}(1),1) scat_locs(sig_grps{grp_ix}(2),1)];
                ys = [scat_locs(sig_grps{grp_ix}(1),2) scat_locs(sig_grps{grp_ix}(2),2)];
                line(xs, ys, 'LineWidth',sig_cnt_roi(grp_ix,roi_ix),'Color',plt.line_color);
                v_stats.ZoneCentroid(grp_ix,:) = [mean(xs) mean(ys)];
            end
        end
        if numel(stat_regs)==3
            v_stats.ZoneCentroid(end,:) = mean(scat_locs,1);
            text_spacers = plt.text_spacers{2};
        else
            text_spacers = plt.text_spacers{1};
        end
        if ~any(sig_cnt_roi(grp_idx,roi_ix))
            text_spacers = zeros(size(text_spacers));
        end
    end
    
    % Plot text for numbers of electrodes in each zone
    for z_ix = 1:n_groups
        text(v_stats.ZoneCentroid(z_ix,1)+text_spacers(z_ix,1),...
             v_stats.ZoneCentroid(z_ix,2)+text_spacers(z_ix,2),...
            num2str(sig_cnt_roi(z_ix,roi_ix)),'FontSize',plt.text_sz,...
            'HorizontalAlignment',plt.text_align);
    end
    axis equal; axis(gca, plt.axis_vis);
    title([atlas_id '-' roi_id '-' roi_list{roi_ix} ': # SBJ = ' num2str(numel(SBJs)) ...
        ', # Sig Elecs = ' num2str(sum(any(sig_roi_all==roi_ix,2))) ' / ' num2str(sum(vertcat(roi_mat{:})==roi_ix)) ...
       ' (' num2str(sum(any(sig_roi_all==roi_ix,2))/sum(vertcat(roi_mat{:})==roi_ix),'%.4f') ')'],...
       'FontSize',plt.title_sz,'Color',roi_colors{roi_ix});
    legend(v,venn_legend,'FontSize',plt.legend_sz);
    
    if save_fig
        fig_fname = [save_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
