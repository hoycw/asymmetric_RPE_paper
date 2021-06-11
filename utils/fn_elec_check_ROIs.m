function fn_elec_check_ROIs(SBJ, reref, varargin)
%% Plot recons for each probe with pial, white matter, and inflated to check ROI assignments

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'proc_id') && ischar(varargin{v+1})
            proc_id = varargin{v+1};
        elseif strcmp(varargin{v},'view_space') && ischar(varargin{v+1})
            view_space = varargin{v+1};
        elseif strcmp(varargin{v},'reg_type') && ischar(varargin{v+1})
            reg_type = varargin{v+1};
        elseif strcmp(varargin{v},'atlas_id') && ischar(varargin{v+1})
            atlas_id = varargin{v+1};
        elseif strcmp(varargin{v},'plot_inflated')
            plot_inflated = varargin{v+1};
        elseif strcmp(varargin{v},'plot_ortho')
            plot_ortho = varargin{v+1};
        elseif strcmp(varargin{v},'probes') && iscell(varargin{v+1})
            probes = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('proc_id','var');       proc_id = 'main_ft'; end
if ~exist('view_space','var');    view_space = 'pat'; end
if ~exist('reg_type','var');      reg_type = ''; end
if ~exist('atlas_id','var');      atlas_id = 'Dx'; end
if ~exist('plot_inflated','var'); plot_inflated = 0; end
if ~exist('plot_ortho','var');    plot_ortho = 1; end

%% Load elec and get ROI labels
% Check which root directory
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end
if reref
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_' atlas_id '_man.mat'];
else
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '.mat'];
end
load(elec_fname);

%% Add fields for orig elecs
if ~reref
    elec.gROI = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'gROI');
    elec.ROI  = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'ROI');
end

% Assign ROIs
elec.color = fn_roi2color(elec.ROI);

%% Load meshes
pial_l = fn_load_recon_mesh(SBJ,view_space,reg_type,'pial','l');
pial_r = fn_load_recon_mesh(SBJ,view_space,reg_type,'pial','r');
if exist([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_wm_lh.mat'],'file')
    wm_l   = fn_load_recon_mesh(SBJ,view_space,reg_type,'wm','l');
    wm_r   = fn_load_recon_mesh(SBJ,view_space,reg_type,'wm','r');
end

%% Plot inflated (no elecs)
if plot_inflated
    % Load inflated surfaces
    infl_l = fn_load_recon_mesh(SBJ,view_space,reg_type,'inflated','l');
    infl_r = fn_load_recon_mesh(SBJ,view_space,reg_type,'inflated','r');
    
    % Plot inflated surfaces
    i_l = figure('Name',[SBJ ' inflated L']);
    ft_plot_mesh(infl_l, 'vertexcolor', 'curv');
    material dull; lighting gouraud;
    l = camlight;
    set(i_l, 'windowkeypressfcn',   @cb_keyboard);
    
    i_r = figure('Name',[SBJ ' inflated L']);
    ft_plot_mesh(infl_r, 'vertexcolor', 'curv');
    material dull; lighting gouraud;
    l = camlight;
    set(i_r, 'windowkeypressfcn',   @cb_keyboard);
end

%% Plot elecs + meshes
% Select probes to plot
if exist('probes','var')
    probe_ix = zeros(size(probes));
    for p = 1:numel(probe_ix)
        probe_ix(p) = find(strcmp(SBJ_vars.ch_lab.probes,probes{p}));
    end
else
    probe_ix = 1:numel(SBJ_vars.ch_lab.probes);
end

% Plot each probe
for p = probe_ix
    cfgs.channel = ft_channelselection([SBJ_vars.ch_lab.probes{p} '*'], elec.label);
    probe = fn_select_elec(cfgs, elec);
    
    % Plot pial
    mesh_alpha = 0.7;
    if numel(probe.label)>1; hemi = probe.hemi{2}; else hemi = probe.hemi{1};end
    pial = figure('Name',[SBJ ' pial ' probe.label{1} ' : ' probe.label{end}],...
                    'Units','normalized');
    set(pial,'OuterPosition',[0 0 0.5 1]);
    ft_plot_mesh(eval(['pial_' hemi]), 'vertexcolor', 'curv', 'facealpha', mesh_alpha);
    for e = 1:numel(probe.label)
        cfgs.channel = probe.label(e);
        tmp = fn_select_elec(cfgs, probe);
        ft_plot_sens(tmp, 'elecshape', 'sphere', 'facecolor', tmp.color, 'label', 'label');
    end
    material dull; lighting gouraud;
    l = camlight;
    set(pial, 'windowkeypressfcn',   @cb_keyboard);
    
    % Plot WM
    if strcmp(SBJ_vars.ch_lab.probe_type{p},'seeg') && exist('wm_l','var')
        wm = figure('Name',[SBJ ' white ' probe.label{1} ' : ' probe.label{end}],...
                    'Units','normalized');
        set(wm, 'OuterPosition', [0.5 0 0.5 1]);
        ft_plot_mesh(eval(['wm_' hemi]), 'vertexcolor', 'curv', 'facealpha', 0.5);
        for e = 1:numel(probe.label)
            cfgs.channel = probe.label(e);
            tmp = fn_select_elec(cfgs, probe);
            ft_plot_sens(tmp, 'elecshape', 'sphere', 'facecolor', tmp.color, 'label', 'label');
        end
        material dull; lighting gouraud;
        l = camlight;
        set(wm, 'windowkeypressfcn',   @cb_keyboard);
    end
end

%% Plot Elecs
if plot_ortho %any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
    if ~reref
        fn_view_recon(SBJ, '', 'ortho', view_space, reg_type, 1, 'b', 1);
    else
        fn_view_recon(SBJ, 'main_ft', 'ortho', view_space, reg_type, 1, 'b', 1);
    end
end

end
