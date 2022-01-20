function SBJ11b_HFA_conn_mLME(SBJs, proc_id, an_id, model_id, conn_id, lme_formula)
%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Data for each subject and stack into tables:
conn_data = [];
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('loading subject %s\n',SBJ)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
    % eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
    %eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
    %if ~strcmp(st.an_style,'mLME'); error('This script is only for mass mLME!'); end

    % Load behavioral data, model, and connectivity data
    load([SBJ_vars.dirs.events,SBJ,'_bhv_',proc_id,'_final.mat']);
    load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
    load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'_', conn_id,'.mat']);
    [reg_lab, ~, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
    
    % select regressors (provisional)
%     model = model(:,1:end-1);
%     reg_lab = reg_lab(1:end-1);
    
    pair_labels = conn.pair_labels;
    try time = conn.time; catch; time = 0; end
    for pl = 1:numel(pair_labels)
        cplab = pair_labels{pl};
        cnbins = size(conn.coeff{pl},2);
        cntimes = size(conn.coeff{pl},5);
        if s == 1
           conn_data{pl} = cell(numel(SBJs), cntimes, cnbins);
        end
        for cb = 1:cnbins
            for cnt = 1:cntimes
                nch1 = numel(conn.chann_labels{pl}{1});
                nch2 = numel(conn.chann_labels{pl}{2});
                conn_data{pl}{s,cnt,cb} = cell(size(conn.coeff{pl},1)*nch1*nch2,6);
                dcount = 0;
                for ch1 = 1:nch1
                    lab1 = conn.chann_labels{pl}{1}{ch1};
                    for ch2 = 1:nch2
                        lab2 = conn.chann_labels{pl}{2}{ch2};
                        cdata = squeeze(conn.coeff{pl}(:,cb,ch1,ch2,cnt));
                        for cdp = 1:numel(cdata)
                            dcount = dcount + 1;
                            conn_data{pl}{s,cnt,cb}{dcount,1} = cdata(cdp);
                            conn_data{pl}{s,cnt,cb}{dcount,2} = SBJ;
                            conn_data{pl}{s,cnt,cb}{dcount,3} = [SBJ '_' lab1 '_to_' lab2];
                            for mm = 1:size(model,2)
                                conn_data{pl}{s,cnt,cb}{dcount,3 + mm} = model(cdp,mm);
                            end
                        end
                    end
                end
            end
        end
    end
    clear dirs_fields field_ix SBJ_id SBJ_vars conn
end

conn_tables = [];
for pl = 1:numel(conn_data)
    conn_tables{pl} = cell(size(conn_data{pl},2),size(conn_data{pl},3));
    for cb = 1:size(conn_data{pl},3)
        for ct = 1:size(conn_data{pl},2)
            ctable = cell2table(vertcat(conn_data{pl}{:,ct,cb}));
            ctable.Properties.VariableNames = [{'y','sub','chan'},reg_lab];
            conn_tables{pl}{ct,cb} = ctable;
        end
    end
end
clear conn_data
%% Run LMEs
LMEs = {};
coefs = {};
lower = {};
upper = {};
pvals = {};
% lme_formula = ['y~sRPE + uRPE + EV + (1 + sRPE + uRPE + EV | sub) +', ...
%                 '(1 + sRPE + uRPE + EV | sub:chan)'];
for pl = 1:numel(conn_tables)
    nb = size(conn_tables{pl},2);
    nt = size(conn_tables{pl},1);
    LMEs{pl} = cell(nt,nb);
    pvals{pl} = NaN(size(model,2) + 1,nt,nb);
    coefs{pl} = NaN(size(model,2) + 1,nt,nb);
    lower{pl} = NaN(size(model,2) + 1,nt,nb);
    upper{pl} = NaN(size(model,2) + 1,nt,nb);
    for b = 1:nb
        for t = 1:nt
            fprintf('fitting label pair %d time point %d and bin %d\n',pl,t,b)
            LMEs{pl}{t,b} = fitlme(conn_tables{pl}{t,b}, lme_formula);
            pvals{pl}(:,t,b) = LMEs{pl}{t,b}.Coefficients.pValue;
            coefs{pl}(:,t,b) = LMEs{pl}{t,b}.Coefficients.Estimate;
            lower{pl}(:,t,b) = LMEs{pl}{t,b}.Coefficients.Lower;
            upper{pl}(:,t,b) = LMEs{pl}{t,b}.Coefficients.Upper;
        end
    end
end
%% FDR correction
qvals = cell(size(pvals));
for pl = 1:numel(pvals)
    qvals{pl} = NaN(size(pvals{pl}));
    for cf = 1:size(pvals{pl},1)
        for b = 1:size(pvals{pl},3)
            [~, ~, ~, qvals{pl}(cf,:,b)] = fdr_bh(squeeze(pvals{pl}(cf,:,b)));
        end
    end
end
% qvals = cell(size(pvals));
% for pl = 1:numel(pvals)
%     qvals{pl} = NaN(size(pvals{pl}));
%     qvals{pl}(1,:,:) = pvals{pl}(1,:,:);
%     for b = 1:size(pvals{pl},3)
%         [~, ~, ~, qvals{pl}(:,:,b)] = fdr_bh(squeeze(pvals{pl}(:,:,b)));
%     end
% end
%% Store in a structure and save
conn_stats =  [];
%conn_stats.LMEs = LMEs;
conn_stats.feature = reg_lab;
conn_stats.coefs = coefs;
conn_stats.lower = lower;
conn_stats.upper = upper;
conn_stats.pvals = pvals;
conn_stats.qvals = qvals;
conn_stats.pair_labels = pair_labels;
conn_stats.time = time;
conn_stats.dimord = 'coef_time_bin';

stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir proc_id '_' model_id '_' an_id '_' conn_id '.mat'];
save(out_fname,'-v7.3','conn_stats')
%% Extract electrode-wise coefficients.
chan_coef = {};
chan_lower = {};
chan_upper = {};
chan_pval = {};
chan_labels = {};
for r = 1:numel(LMEs)
    chan_coef{r} = NaN(size(model,2) + 1, length(unique(conn_tables{r}{1,1}.chan)),...
                    size(LMEs{r},1),size(LMEs{r},1));
    chan_lower{r} = chan_coef{r};
    chan_upper{r} = chan_coef{r};
    chan_pval{r} = chan_coef{r};
    for t = 1:size(LMEs{r},1)
        for b = 1:size(LMEs{r},2)
            [~,~,rndstats] =  randomEffects(LMEs{r}{t,b});
            
            %select channel random effect
            rndstats = rndstats(strcmp(rndstats.Group,'sub:chan'),:);
            if t == 1 & b == 1
                chan_labels{r} = unique(rndstats.chan); 
            end
            chan_coef{r}(:,:,t,b) = reshape(rndstats.Estimate,size(chan_coef{r},1:2));
            chan_lower{r}(:,:,t,b) = reshape(rndstats.Lower,size(chan_coef{r},1:2));
            chan_upper{r}(:,:,t,b) = reshape(rndstats.Upper,size(chan_coef{r},1:2));
            chan_pval{r}(:,:,t,b) = reshape(rndstats.pValue,size(chan_coef{r},1:2));
        end
    end
    chan_coef{r} = permute(chan_coef{r},[2,1,3,4]);
    chan_lower{r} = permute(chan_lower{r},[2,1,3,4]);
    chan_upper{r} = permute(chan_upper{r},[2,1,3,4]);
    chan_pval{r} = permute(chan_pval{r},[2,1,3,4]);
end

% FDR correction
chan_qval = {};
for r = 1:numel(chan_pval)
    for ch = 1:size(chan_pval{r},1)
        for rg = 1:size(chan_pval{r},2)
            for b = 1:size(chan_pval,4)
                [~,~,~,chan_qval{r}(ch,rg,:,b)] = fdr_bh(chan_pval{r}(ch,rg,:,b));
            end
        end
    end
end

% Store in a structure
conn_stats_chan =  [];
%conn_stats.LMEs = LMEs;
conn_stats_chan.feature   = reg_lab;
conn_stats_chan.time      = st.win_center;
conn_stats_chan.coefs = chan_coef;
conn_stats_chan.lower = chan_lower;
conn_stats_chan.upper = chan_upper;
conn_stats_chan.pvals = chan_pval;
conn_stats_chan.qvals = chan_qval;
conn_stats_chan.label = roi_list;
conn_stats_chan.chan_label = chan_labels;
conn_stats_chan.win_lim   = win_lim;
conn_stats_chan.win_lim_s = hfa_time(win_lim);
conn_stats_chan.dimord = 'chan_coef_time_bin';

% save
fprintf('=================== Saving channel effects ======================\n');
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir model_id '_' stat_id '_' an_id '_conn_chancoef.mat'];
save(out_fname,'-v7.3','conn_stats_chan')
