function [design,col_names] = fn_build_design_predictors(SBJ,stat_id)
%% Build design matrix of predictors for regression analysis
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Variable Handline
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

%% Load Data
load([SBJ_vars.dirs.events SBJ '_trl_info_final.mat']);
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

%% Build Design Matrix
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
groups = cell([1 numel(grp_lab)]);
levels = cell([1 numel(grp_lab)]);
for grp_ix = 1:numel(grp_lab)
    [levels{grp_ix}, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
    cond_idx = fn_condition_index(grp_lab{grp_ix}, trl_info);
    groups{grp_ix} = cell(size(cond_idx));
    for level_ix = 1:numel(levels{grp_ix})
        groups{grp_ix}(cond_idx==level_ix) = levels{grp_ix}(level_ix);
    end
end

% Add hit/miss
hit_str = cell(size(trl_info.hit));
hit_str(logical(trl_info.hit))  = {'hit'};
hit_str(~logical(trl_info.hit)) = {'miss'};

%% Create matrix with trl_info fields
col_names = {'trl_n','blk','rt','tol'};%'hit' being handled as categorical string
design = zeros([numel(trl_info.trl_n) numel(col_names)]);
for c_ix = 1:numel(col_names)
    design(:,c_ix) = trl_info.(col_names{c_ix});
end

%% Add signed quantitative surprise (inverse distance from tolerance)
[tim_lab,~,~,~] = fn_condition_label_styles('Tim');
tim_idx = fn_condition_index('Tim', trl_info);

% Compute distance from closest tolearnce bound
%   Positive = win, Negative = loss
signed_pe = zeros(size(trl_info.rt));
for cond_ix = 1:numel(tim_lab)
    if strcmp(tim_lab{cond_ix},'Er')
        signed_pe(tim_idx==cond_ix) = trl_info.rt(tim_idx==cond_ix)-(prdm_vars.target-trl_info.tol(tim_idx==cond_ix));
    else
        signed_pe(tim_idx==cond_ix) = (prdm_vars.target+trl_info.tol(tim_idx==cond_ix))-trl_info.rt(tim_idx==cond_ix);
    end
end

% Invert to get prediction error
signed_pe = signed_pe.^-1;
col_names{end+1} = 'signPE';
design = horzcat(design,signed_pe);

% Convert to pure prediction error
col_names{end+1} = 'absPE';
design = horzcat(design,abs(signed_pe));

%% Add previous trial fields
% Previous trial outcome + surprise
out_ix   = find(strcmp(grp_lab,'Out'));
sur_ix   = find(strcmp(grp_lab,'Dif*Out'));
prev_out = cell(size(groups{out_ix}));
prev_sur = cell(size(groups{sur_ix}));
prev_out{1} = 'NA';
prev_sur{1} = 'NA';
for t_ix = 2:numel(prev_out)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_out(t_ix) = groups{out_ix}(t_ix-1);
        prev_sur(t_ix) = groups{sur_ix}(t_ix-1);
    else
        prev_out{t_ix} = 'NA';
        prev_sur{t_ix} = 'NA';
    end
end

% Previous RT
prev_rt = nan(size(trl_info.rt));
for t_ix = 2:numel(prev_rt)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_rt(t_ix) = trl_info.rt(t_ix-1);
    end
end

% Previous Tolerance
prev_tol = nan(size(trl_info.tol));
for t_ix = 2:numel(prev_tol)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_tol(t_ix) = trl_info.tol(t_ix-1);
    end
end

% Previous trial PE
prev_spe = nan(size(signed_pe));
for t_ix = 2:numel(prev_spe)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_spe(t_ix) = signed_pe(t_ix-1);
    end
end

% Add these to predictor matrix
col_names{end+1} = 'n1_rt';
design = horzcat(design,prev_rt);
col_names{end+1} = 'n1_tol';
design = horzcat(design,prev_tol);
col_names{end+1} = 'n1_signPE';
design = horzcat(design,prev_spe);
col_names{end+1} = 'n1_absPE';
design = horzcat(design,abs(prev_spe));

%% Combine categorical and continuous variables
% Convert continuous to cell
design = num2cell(design);

% Combine column names
col_names = [grp_lab {'hit', 'n1_Out', 'n1_Dif*Out'} col_names];
design = horzcat(groups{:},hit_str,prev_out,prev_sur,design);

% Convert NaN to 'NA'
for r = 1:size(design,1)
    for c = 1:size(design,2)
        if any(isnan(design{r,c})); design{r,c} = 'NA'; end
    end
end

end