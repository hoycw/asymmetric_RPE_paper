function [pred,col_names,col_vals] = fn_build_design_predictors(SBJ,stat_id)
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
design = cell([1 numel(grp_lab)]);
levels = cell([1 numel(grp_lab)]);
for grp_ix = 1:numel(grp_lab)
    [levels{grp_ix}, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
    design{grp_ix} = fn_condition_index(grp_lab{grp_ix}, trl_info);
end

%% Create matrix with trl_info fields
col_names = {'trl_n','blk','hit','rt','tol'};
col_vals  = cell(size(col_names));
pred = zeros([numel(trl_info.trl_n) numel(col_names)]);
for c_ix = 1:numel(col_names)
    if strcmp(col_names{c_ix},'cond')
        pred(:,c_ix) = cond_idx;
    else
        pred(:,c_ix) = trl_info.(col_names{c_ix});
    end
    switch col_names{c_ix}
        case {'trl_n','blk'}
            col_vals{c_ix} = 'category';
        case {'rt','tol'}
            col_vals{c_ix} = 'continuous';
        case 'cond'
            col_vals{c_ix} = cond_lab;
        case 'hit'
            col_vals{c_ix} = {'Ls','Wn'};
        otherwise
            error(['unknown column: ' col_names{c_ix}]);
    end
end

%% Add in design
col_names = [grp_lab col_names];
col_vals = [levels col_vals];
pred = horzcat(design{:},pred);

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
col_names{end+1} = 'sPE';
col_vals{end+1}  = 'continuous';
pred = horzcat(pred,signed_pe);

% Convert to pure prediction error
col_names{end+1} = 'uPE';
col_vals{end+1}  = 'continuous';
pred = horzcat(pred,abs(signed_pe));


%% Add previous trial fields
% Previous RT
prev_rt = nan(size(trl_info.rt));
for t_ix = 2:numel(prev_rt)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_rt(t_ix) = trl_info.rt(t_ix-1);
    end
end

col_names{end+1} = 'rtn-1';
col_vals{end+1} = 'continuous';
pred = horzcat(pred,prev_rt);

% Previous trial outcome + surprise
out_ix   = find(strcmp(grp_lab,'Out'));
sur_ix   = find(strcmp(grp_lab,'Dif*Out'));
prev_out = nan(size(design{out_ix}));
prev_sur = nan(size(design{sur_ix}));
for t_ix = 2:numel(prev_out)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_out(t_ix) = design{out_ix}(t_ix-1);
        prev_sur(t_ix) = design{sur_ix}(t_ix-1);
    end
end

col_names{end+1} = 'Outn-1';
col_vals{end+1} = levels{out_ix};
pred = horzcat(pred,prev_out);
col_names{end+1} = 'D*On-1';
col_vals{end+1} = levels{sur_ix};
pred = horzcat(pred,prev_sur);

% Previous trial PE
prev_spe = nan(size(signed_pe));
for t_ix = 2:numel(prev_spe)
    if trl_info.blk(t_ix)==trl_info.blk(t_ix-1) && trl_info.trl_n(t_ix)==trl_info.trl_n(t_ix-1)+1
        prev_spe(t_ix) = signed_pe(t_ix-1);
    end
end

col_names{end+1} = 'sPEn-1';
col_vals{end+1} = 'continuous';
pred = horzcat(pred,prev_spe);
col_names{end+1} = 'uPEn-1';
col_vals{end+1} = 'continuous';
pred = horzcat(pred,abs(prev_spe));


end