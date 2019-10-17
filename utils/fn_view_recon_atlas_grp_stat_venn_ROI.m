function fn_view_recon_atlas_grp_stat_venn_ROI(SBJs, proc_id, stat_id, an_id, reg_type, show_labels,...
                            hemi, atlas_id, roi_id, plot_roi, mirror, varargin)
%% Plot a reconstruction with electrodes colored according to statistics
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pcon': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Process plotting params
% Error cases
if strcmp(hemi,'b') && ~strcmp(plot_roi,'OFC')
    error('hemi must be l or r for all non-OFC plots');
end
if ~any(strcmp(plot_roi,{'LPFC','MPFC','INS','OFC','TMP','PAR','lat','deep'}))
    error('roi_id needs to be a lobe, "lat", or "deep"');
end

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype')
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Implement the default options
% Add default view_angle if not defined
if ~exist('view_angle','var')
    view_angle = fn_get_view_angle(hemi,plot_roi);
    view_str = 'def';
end
% Adjust view angle if custom
if ischar(view_angle)
    view_str = view_angle;
    if (strcmp(hemi,'l') && strcmp(view_angle,'med')) || (strcmp(hemi,'r') && strcmp(view_angle,'lat'))
        view_angle = [90 0];
    elseif (strcmp(hemi,'l') && strcmp(view_angle,'lat')) || (strcmp(hemi,'r') && strcmp(view_angle,'med'))
        view_angle = [-90 0];
    end
end
if ~exist('mesh_alpha','var')
    mesh_alpha = 0.3;   % assume SEEG
end
if mirror && strcmp(hemi,'b')
    error('why mirror if hemi b?');
end
if ~exist('save_fig','var')
    save_fig = 0;
end
if ~exist('fig_ftype','var')
    fig_ftype = 'fig';
end
if show_labels
    lab_arg = 'label';
else
    lab_arg = 'off';
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

% ROI info
if any(strcmp(plot_roi,{'deep','lat'}))
    [plot_roi_list, ~] = fn_roi_label_styles(plot_roi);
else
    plot_roi_list = {plot_roi};
end

if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

%% Load data
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
if numel(grp_lab) < 2 || numel(grp_lab) > 3; error('why venn?'); end
venn_colors = fn_venn_colors(numel(grp_lab));
all_color   = [1 1 1];

elec_sbj    = cell([numel(SBJs) 1]);
elec_sig    = cell([numel(SBJs) 1]);
good_sbj    = true([numel(SBJs) 1]);
sig_roi_mat = cell([numel(SBJs) 1]);
roi_mat     = cell([numel(SBJs) 1]);
all_roi_colors = cell([numel(grp_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    
    % Append SBJ name to labels
    orig_labels = elec_sbj{sbj_ix}.label;
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
    end
    
    % Select ROI mesh matches
    plot_elecs = zeros([numel(elec_sbj{sbj_ix}.label) numel(plot_roi_list)]);
    for roi_ix = 1:numel(plot_roi_list)
        plot_elecs(:,roi_ix) = strcmp(elec_sbj{sbj_ix}.(roi_field),plot_roi_list{roi_ix});
    end
    % Remove electrodes that aren't in atlas ROIs & hemisphere
    if mirror
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, 'b', atlas_id, roi_id);
        hemi_str = [hemi 'b'];
    else
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
        hemi_str = hemi;
    end
    plot_roi_elecs = intersect(roi_elecs, elec_sbj{sbj_ix}.label(any(plot_elecs,2)));
    
    % Match elecs to plotting ROIs
    roi_mat{sbj_ix} = zeros([numel(elec_sbj{sbj_ix}.label) 1]);
    for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
        if any(strcmp(elec_sbj{sbj_ix}.label{ch_ix},plot_roi_elecs))
            roi_mat{sbj_ix}(ch_ix) = find(strcmp(plot_roi_list,elec_sbj{sbj_ix}.(roi_field){ch_ix}));
        end
    end
    
    %% Load Stats
    sig_mat = zeros([numel(elec_sbj{sbj_ix}.label) numel(grp_lab)]);
    if any(roi_mat{sbj_ix})
        % Mirror hemispheres
        if mirror
            elec_sbj{sbj_ix}.chanpos(~strcmp(elec_sbj{sbj_ix}.hemi,hemi),1) = ...
                -elec_sbj{sbj_ix}.chanpos(~strcmp(elec_sbj{sbj_ix}.hemi,hemi),1);
        end
        
        load([SBJ_vars.dirs.proc SBJ '_nANOVA_ROI_' stat_id '_' an_id '.mat'],'w2');
        % HACK!!! remove some to match with PRJ_Stroop elec files
        cfgs = [];
        if strcmp(SBJ,'CP24')
            cfgs.channel = {'all','-RTO4'};
        elseif strcmp(SBJ,'IR57')
            cfgs.channel = {'all','-LAM4-5','-LAM5-6','-RIN8-9','-RIN9-10',...
                '-RSM1-2','-RSM2-3','-RSM3-4','-RTI9-10'};
        elseif strcmp(SBJ,'IR68')
            cfgs.channel = {'all','-LPC5-6','-LPC6-7','-LPC7-8'};
        end
        w2 = ft_selectdata(cfgs,w2);
        if numel(elec_sbj{sbj_ix}.label)~=numel(w2.label) || ~all(strcmp(orig_labels,w2.label))
            error('mismatched labels in elec and w2!');
        end
        for ch_ix = 1:numel(w2.label)
            % FDR correct pvalues for ANOVA
            [~, ~, ~, qvals] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));
            % Consolidate to binary sig/non-sig
            for grp_ix = 1:numel(grp_lab)
                if any(qvals(grp_ix,:)<=0.05,2)
                    sig_mat(ch_ix,grp_ix) = 1;
                end
            end
        end
        
        %% Compile Statistics (if any significance in plot_roi)
        if any(any(sig_mat(roi_mat{sbj_ix}~=0,:),2))
            % Print # and % sig
            print_nums = zeros([2 numel(grp_lab)]);
            for st_ix = 1:numel(grp_lab)
                print_nums(1,st_ix) = sum(sig_mat(:,st_ix));
                print_nums(2,st_ix) = sum(sig_mat(:,st_ix))/size(sig_mat,1);
            end
            fprintf(['%-10s' repmat('%-3i(%.3f)\t',[1 numel(grp_lab)]) '\n'],[SBJ ' sig:'],print_nums(:));
            
            % Keep intersection of significant and ROI matched electrodes
            sig_idx = any(sig_mat,2);
            roi_idx = roi_mat{sbj_ix}~=0;
            sig_roi_mat{sbj_ix} = sig_mat(roi_idx & sig_idx,:).*roi_mat{sbj_ix}(roi_idx & sig_idx);
            % Select sig elecs for plotting
            cfgs = []; cfgs.channel = elec_sbj{sbj_ix}.label(sig_idx & roi_idx);
            elec_sig{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});
            
            % For significant elecs in plot_roi, select color and add to elec
            all_roi_colors{sbj_ix} = zeros([numel(elec_sig{sbj_ix}.label) 3]);
            for sig_ix = 1:numel(elec_sig{sbj_ix}.label)
                grp_ix = find(sig_roi_mat{sbj_ix}(sig_ix,:));
                if numel(grp_ix)==numel(grp_lab) && numel(grp_lab) == 3
                    all_roi_colors{sbj_ix}(sig_ix,:) = all_color;
                elseif numel(grp_ix)==2
                    all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{grp_ix(1),grp_ix(2)};
                else
                    all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{grp_ix,grp_ix};
                end
            end
            
            fprintf('\t%s has %i sig channels in %s hemi %s\n',SBJ,size(sig_roi_mat{sbj_ix},1),atlas_id,hemi);
        else
            % Print no significant elecs
            good_sbj(sbj_ix) = false;
            fprintf(2,'\t%s has %i channels in %s hemi %s, but none are significant\n',...
                SBJ,sum(roi_mat{sbj_ix}~=0),atlas_id,hemi);
        end
    else
        % Print no ROI match
        good_sbj(sbj_ix) = false;
        fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJ,atlas_id,hemi_str);
    end
    clear SBJ SBJ_vars SBJ_vars_cmd w2
end

%% Combine elec structs
elec = ft_appendsens([],elec_sig{good_sbj});
elec.color = vertcat(all_roi_colors{:});    % appendsens strips that field

%% Load Atlas
atlas = fn_load_recon_atlas([],atlas_id);

% Get Atlas-ROI mapping
atlas_labels = fn_atlas_roi_select_mesh(atlas_id, plot_roi, hemi);

% Can't plot unconnected meshes (I think), so create two meshes
if strcmp(plot_roi,'OFC')
    % Treat R hemi as new ROI
    r_ix = ~cellfun(@isempty,strfind(atlas_labels,'rh'));
    r_labels = atlas_labels(r_ix);
    atlas_labels = atlas_labels(~r_ix);
elseif strcmp(plot_roi,'deep')
    mtl_ix = ~cellfun(@isempty,strfind(atlas_labels,'Hippocampus')) | ...
             ~cellfun(@isempty,strfind(atlas_labels,'Amygdala'));
%     mtl_labels = atlas_labels(mtl_ix);
    atlas_labels = atlas_labels(~mtl_ix);
end

%% Select ROI mesh
cfg = [];
cfg.inputcoord = atlas.coordsys;
cfg.atlas = atlas;
cfg.roi = atlas_labels;
roi_mask = ft_volumelookup(cfg,atlas);
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = roi_mask;

if exist('r_labels','var')
    cfg.roi = r_labels;
    r_mask  = ft_volumelookup(cfg,atlas);
    
    r_seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    r_seg.brain = r_mask;
elseif exist('mtl_labels','var')
    cfg.roi  = mtl_labels;
    mtl_mask = ft_volumelookup(cfg,atlas);
    
    mtl_seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    mtl_seg.brain = mtl_mask;
end

cfg = [];
cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 100000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
roi_mesh = ft_prepare_mesh(cfg, seg);
if exist('r_seg','var')
    r_mesh = ft_prepare_mesh(cfg, r_seg);
elseif exist('mtl_seg','var')
    mtl_mesh = ft_prepare_mesh(cfg, mtl_seg);
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
if save_fig
    out_dir = [root_dir 'PRJ_Error/results/HFA/GRP_recon_venn/' stat_id '/' an_id '/'];
    if ~exist(out_dir,'dir')
        [~] = mkdir(out_dir);
    end
end
plot_name = ['GRP_' strjoin(grp_lab,'-') '_' stat_id '_' an_id '_' plot_roi '_' hemi_str '_' view_str];
f = figure('Name',plot_name);

% Plot 3D mesh
ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
if exist('r_mesh','var')
    ft_plot_mesh(r_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
elseif exist('mtl_mesh','var')
    ft_plot_mesh(mtl_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
end

% Plot electrodes on top
for e = 1:numel(elec{grp_ix}.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(f(grp_ix), 'windowkeypressfcn',   @cb_keyboard);

if save_fig
    fig_fname = [out_dir plot_name fig_ftype];
    fig_fname = strrep(fig_fname,'*','x');
    saveas(f(grp_ix),fig_fname);
end

%% Plot legend in separate figure
venn_name = [plot_name '_leg'];
v = figure('Name',venn_name); hold on;
venn_legend = {}; leg_ix = 0;
for grp_ix1 = 1:numel(grp_ix)
    for grp_ix2 = grp_ix1:numel(grp_lab)
        scatter(grp_ix1,grp_ix2,500,'MarkerFaceColor',venn_colors{grp_ix1,grp_ix2},...
                'MarkerEdgeColor','k');
        leg_ix = leg_ix + 1;
        if grp_ix1==grp_ix2
            venn_legend{leg_ix} = [grp_lab{grp_ix1} num2str(grp_ix1)];
        else
            venn_legend{leg_ix} = [grp_lab{grp_ix1} num2str(grp_ix1) '+'...
                                   grp_lab{grp_ix2} num2str(grp_ix2)];
        end
    end
end
if numel(grp_lab)>2
    scatter(numel(grp_lab),1,500,'MarkerFaceColor','w','MarkerEdgeColor','k');
    venn_legend{leg_ix+1} = strjoin(grp_lab,'+');
end
xlim([0 numel(grp_lab)+1]);
ylim([0 numel(grp_lab)+1]);
legend(venn_legend);

if save_fig
    saveas(v, [out_dir venn_name '.' fig_ftype]);
end

end
