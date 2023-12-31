function SBJ08g_HFA_crRT_mLME(SBJs, proc_id, an_id, model_id, stat_id, atlas_id, roi_id, lme_formula)
%% function SBJ08a_crRT_mGLM(SBJ,an_id,stat_id)
%   Run GLM with sliding windows for given time-frequency analysis
%   Correlation with RT and regression of RT as confound supported
%   Model based on log fit to behavior:
%       reward prediction, signed PE, unsigned PE
% INPUTS:
%   SBJ [str] - subject ID
%   proc_id [str] - ID of processing pipeline
%   an_id [str] - ID of HFA analysis to run stats
%   model_id [str] - ID of the model used in GLM
%   stat_id [str] - ID of the statistical parameters to extract from HFA
% OUTPUTS:
%   beta [struct] - pseudo-FT structure with main GLM outputs
%   rt [FT struct] - output of correlation with RT if st.rt_corr==1
%   st [struct] - stat params loaded via stat_id

%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];

% rng('shuffle'); % seed randi with time

%% Load Data
hfa_data = [];
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('loading subject %s\n',SBJ)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
    % eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
    if ~strcmp(st.an_style,'mLME'); error('This script is only for mass mLME!'); end
    
    % Load behavioral data, model, and HFA
    load([SBJ_vars.dirs.events,SBJ,'_bhv_',proc_id,'_final.mat']);
    load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
    load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'.mat']);
    
    % Check if more than one frequency, error for now
    if numel(hfa.freq)>1
        error('HFA has more than one frequency, can''t run on that for now!');
    end
    %% Select channels in ROIs
    [roi_list,~, roi_field] = fn_roi_label_styles(roi_id);
    load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    r_ix = [];
    for r = 1:numel(roi_list)
        cr_ix = find(strcmp(elec.(roi_field), roi_list{r}));
        r_ix = [r_ix; cr_ix];
    end
    %% Select conditions of interest
    [reg_lab, ~, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
    [cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);
    [mdl_cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(mdl.model_cond);
    
    if ~all(contains(cond_lab,mdl_cond_lab))
        error([model_id ' does not cover all conditions in ' stat_id '!']);
    end
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    bhv = fn_select_bhv(bhv, full_cond_idx);
    model = model(full_cond_idx~=0,:);
    
    % select regressors (provisional)
    %     model = model(:,1:end-1);
    %     reg_lab = reg_lab(1:end-1);
    %     disp(reg_lab)
    
    % Select data in stat window
    cfg_trim = [];
    cfg_trim.latency = st.stat_lim;
    cfg_trim.trials  = full_cond_idx~=0;
    cfg_trim.channel =  elec.label(r_ix);
    hfa = ft_selectdata(cfg_trim,hfa);
    hfa_time = hfa.time;
    
    roi_labs = elec.(roi_field)(r_ix);
    chan_labs = elec.label(r_ix);
    %% Compute max z score per condition (for later exclusions)
    % cond_idx = fn_condition_index(cond_lab, bhv);
    %
    % max_z = zeros(size(hfa.label));
    % cfg_avg = [];
    % cfg_avg.avgoverrpt = 'yes';
    % for cond_ix = 1:numel(cond_lab)
    %     % Average within condition
    %     cfg_avg.trials = find(cond_idx==cond_ix);
    %     tmp_avg = ft_selectdata(cfg_avg,hfa);
    %     % Binarize channels by average HFA vs. threshold
    %     max_z = max(max_z,max(squeeze(tmp_avg.powspctrm),[],2));
    % end
    
    %% Run correlations with RT
    if st.rt_corr
        fprintf('================== Running Correlation with RT =======================\n');
        cfg_rt.design           = zscore(bhv.rt);
        cfg_rt.ivar             = 1;
        rt = ft_freqstatistics(cfg_rt, hfa);
    end
    %% Regress off Reaction Time
    if st.regress_rt
        fprintf('================== Regressing RT =======================\n');
        cfg_conf = [];
        cfg_conf.output = 'residual';
        cfg_conf.confound = zscore(bhv.rt);
        hfa = ft_regressconfound(cfg_conf, hfa);
    end
    %% Average HFA in Sliding Windows
    fprintf('================== Averaging HFA within Windows =======================\n');
    win_lim = fn_get_win_lim_from_center(hfa.time,st.win_center,st.win_len);
    
    % Average in windows
    hfa_win = zeros([size(hfa.powspctrm,1) size(hfa.powspctrm,2) size(win_lim,1)]);%size(hfa.powspctrm,3)
    for w_ix = 1:size(win_lim,1)
        hfa_win(:,:,w_ix) = squeeze(nanmean(hfa.powspctrm(:,:,1,win_lim(w_ix,1):win_lim(w_ix,2)),4));
    end
    %% Store data in long format per roi
    fprintf('==================== Converting to long format ========================\n');
    ntimes = size(hfa_win,3);
    
    for r = 1:numel(roi_list)
        croi = roi_list{r};
        if s == 1
            hfa_data{r} = cell(numel(SBJs), ntimes);
        end
        cur_chan = chan_labs(strcmp(roi_labs,croi));
        nch = numel(cur_chan);
        for ct = 1:ntimes
            hfa_data{r}{s,ct} = cell(size(hfa_win,1)*nch, 3 + size(model,2));
            dcount = 0;
            for ch = 1:nch
                lab = cur_chan{ch};
                cdata = squeeze(hfa_win(:,ch,ct));
                for cdp = 1:numel(cdata)
                    dcount = dcount + 1;
                    hfa_data{r}{s,ct}{dcount,1} = cdata(cdp);
                    hfa_data{r}{s,ct}{dcount,2} = SBJ;
                    hfa_data{r}{s,ct}{dcount,3} = lab;
                    for mm = 1:size(model,2)
                        hfa_data{r}{s,ct}{dcount, 3 + mm} = model(cdp,mm);
                    end
                end
            end
        end
    end
    clear dirs_fields field_ix SBJ_id SBJ_vars hfa hfa_win
end

%% put together in a table
hfa_tables = [];
for r = 1:numel(hfa_data)
    hfa_tables{r} = cell(size(hfa_data{r},2),1);
    for ct = 1:size(hfa_data{r},2)
        ctable = cell2table(vertcat(hfa_data{r}{:,ct}));
        ctable.Properties.VariableNames = [{'y','sub','chan'},reg_lab];
        if st.rank_transform == 1
            disp('transforming to rankit')
            ctable.y = norminv((tiedrank(ctable.y)-.5) / height(ctable))*std(ctable.y)+mean(ctable.y);
        end
        % mean center
%         for rl = 1:length(reg_lab)
%              creg = reg_lab{rl};
%              ctable.(creg) = ctable.(creg) - mean(ctable.(creg));
%         end
        hfa_tables{r}{ct} = ctable;
    end
end
clear hfa_data

%% Run LMEs
fprintf('==================== Fitting LME ========================\n');
options = statset('LinearMixedModel');
if st.robust == 1 
    disp('using robust estimation')
    options.RobustWgtFun = 'bisquare';
end

LMEs = {};
coefs = {};
lower = {};
upper = {};
pvals = {};
AIC = {};
BIC = {};
MSE = {};
r2 = {};
resid = {};
fittedv = {};

for r = 1:numel(hfa_tables)
    nt = size(hfa_tables{r},1);
    LMEs{r} = cell(nt,1);
    pvals{r} = NaN(size(model,2) + 1,nt);
    coefs{r} = NaN(size(model,2) + 1,nt);
    lower{r} = NaN(size(model,2) + 1,nt);
    upper{r} = NaN(size(model,2) + 1,nt);
    AIC{r} = NaN(1,nt);
    BIC{r} = NaN(1,nt);
    LL{r} = NaN(1,nt);
    MSE{r} = NaN(1,nt);
    r2{r} = NaN(1,nt);
    resid{r} = NaN(height(hfa_tables{r}{1}),nt);
    fittedv{r} = NaN(height(hfa_tables{r}{1}),nt);
    
    for t = 1:nt
        
        fprintf('fitting label %d and time point %d\n',r,t)
        
        ctable = hfa_tables{r}{t};
        LMEs{r}{t} = fitlme(ctable, lme_formula,'OptimizerOptions',options, 'FitMethod','REML');
        pvals{r}(:,t) = LMEs{r}{t}.Coefficients.pValue;
        coefs{r}(:,t) = LMEs{r}{t}.Coefficients.Estimate;
        lower{r}(:,t) = LMEs{r}{t}.Coefficients.Lower;
        upper{r}(:,t) = LMEs{r}{t}.Coefficients.Upper;
        AIC{r}(1,t) = LMEs{r}{t}.ModelCriterion.AIC;
        BIC{r}(1,t) = LMEs{r}{t}.ModelCriterion.BIC;
        LL{r}(1,t) = LMEs{r}{t}.LogLikelihood;
        MSE{r}(1,t) = LMEs{r}{t}.MSE;
        r2{r}(1,t) = LMEs{r}{t}.Rsquared.Adjusted;
        resid{r}(:,t) = residuals(LMEs{r}{t});
        fittedv{r}(:,t) = fitted(LMEs{r}{t});
    end
end
%% FDR correction
fprintf('==================== FDR correction ========================\n');
qvals = cell(size(pvals));
for r = 1:numel(pvals)
    qvals{r} = NaN(size(pvals{r}));
    for cf = 1:size(pvals{r},1)
        [~, ~, ~, qvals{r}(cf,:)] = fdr_bh(squeeze(pvals{r}(cf,:)));
    end
end
% for r = 1:numel(qvals)
%     qvals{r} = NaN(size(pvals{r}));
%     qvals{r}(1,:) = pvals{r}(1,:); % don't include intercept in correction
%     [~, ~, ~, qvals{r}(2:end,:)] = fdr_bh(pvals{r}(2:end,:));
% end
%% Store in a structure
beta =  [];
%conn_stats.LMEs = LMEs;
beta.feature   = reg_lab;
beta.time      = st.win_center;
beta.coefs = coefs;
beta.lower = lower;
beta.upper = upper;
beta.pvals = pvals;
beta.qvals = qvals;
beta.AIC = AIC;
beta.BIC = BIC;
beta.LL = LL;
beta.MSE = MSE;
beta.r2 = r2;
beta.label = roi_list;
beta.win_lim   = win_lim;
beta.win_lim_s = hfa_time(win_lim);
beta.dimord = 'coef_time';

%% Plot residuals qqplot
for r = 1:numel(hfa_tables)
    [~,lpval,~,~] = lillietest(reshape(resid{r},[numel(resid{r}),1]));

    cf = figure('units','normalized','outerposition',[0 0 1 1],...
                'PaperOrientation','Landscape');
    subplot(1,2,1)
    qqplot(reshape(resid{r},[numel(resid{r}),1])./std(resid{r},[],'all'))
    %ylim([-5,5])
    title(sprintf('Lillefors pval = %.3f',lpval))

    subplot(1,2,2)
    [f, x] = hist(resid{r}, 50); % Create histogram from a normal distribution.
    g = 1 / sqrt(2 * pi) * exp(-0.5 * (x/std(resid{r},[],'all') ) .^ 2); % pdf of the normal distribution

    bar(x, f / trapz(x, f)); hold on
    plot(x, g / trapz(x, g), 'r'); hold off  
    xlim([-5*std(resid{r},[],'all'), 5*std(resid{r},[],'all')])

    %histfit(reshape(resid{r},[numel(resid{r}),1]))
    title(sprintf('median = %.3f',median(resid{r},'all')))
    
    plot_fname1 = [fig_dir proc_id '_' model_id '_' an_id '_' stat_id '_hfa_'...
                  beta.label{r} '_residuals_qqplot.pdf'];
    
    print(plot_fname1,cf,'-dpdf','-fillpage')
end
%close all
%% Residuals vs fitted
for r = 1:numel(hfa_tables)
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
                'PaperOrientation','Landscape');
    s = scatter(reshape(fittedv{r},[numel(fittedv{r}),1]),...
            reshape(resid{r},[numel(resid{r}),1]),'.k');
    s.MarkerEdgeAlpha = .05;
    s.MarkerFaceAlpha = .05;
    title(beta.label{r})
    xlabel('fitted')
    ylabel('residuals')
    plot_fname2 = [fig_dir proc_id '_' model_id '_' an_id '_' stat_id '_hfa_'...
                  beta.label{r} '_residuals_vs_fitted.png'];
    
    print(plot_fname2,cf,'-dpng')%,'-fillpage')
end
close all
%% Residuals vs predictors
for r = 1:numel(hfa_tables)
    cf = figure('units','normalized','outerposition',[0 0 1 1],...
                'PaperOrientation','Landscape');
    for creg = 1:length(reg_lab) 
        subplot(1,length(reg_lab),creg)
        s = scatter(repmat(hfa_tables{r}{1}.(reg_lab{creg}),[size(resid{r},2),1]),...
             reshape(resid{r},[numel(resid{r}),1]),'.k');
        xlabel(reg_lab{creg})
        ylabel('residuals')
        title(reg_lab{creg})
        s.MarkerEdgeAlpha = .05;
        s.MarkerFaceAlpha = .05;
    end
    
    plot_fname3 = [fig_dir proc_id '_' model_id '_' an_id '_' stat_id '_hfa_'...
                  beta.label{r} '_residuals_vs_predictors.png'];
    
    print(plot_fname3,cf,'-dpng')%,'-fillpage')
end
%%
close all
%% Save
fprintf('==================== Saving results ========================\n');
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa.mat'];
%out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_REML.mat'];
save(out_fname,'-v7.3','beta')

%% Extract electrode-wise coefficients.
chan_coef = {};
% chan_lower = {};
% chan_upper = {};
chan_pval = {};
chan_labels = {};

for r = 1:numel(LMEs)
    for t = 1:numel(LMEs{r})
        [~,~,re] = randomEffects(LMEs{r}{t},'DFMethod','residual');%'satterthwaite')%;
        [~,~,fe] = fixedEffects(LMEs{r}{t});
        fprintf('estimating random effects for region %d and timepoint %d\n',r,t)
        
        nr = length(re.Estimate);
        nf = length(fe.Estimate);

        %Compute additive estimates and p-values
        comp_Estimate = nan(nr,1);
        comp_pValue = nan(nr,1);
%        rlevels = unique(rnames.Level);
%             for f = 1:height(fnames)
%                 cfname = fnames.Name{f};
%                 disp(cfname)
%            fmat = zeros(1,height(fnames));
%                 fmat(f) = 1;
        cfnames=fe.Name;
        crnames=re.Name;
        crlevels=re.Level;
        for l = 1:nr
            crname = crnames{l};
            fmat = zeros(1,nf);
            rmat = zeros(1,nr);
            fmat(1,strcmp(cfnames,crname)) = 1;
            rmat(1,l) = 1;
            clevel = crlevels{l};
            %disp(clevel)
%            cslevels = split(clevel,' ');
%             if length(cslevels) > 1
%                 rmat(1,strcmp(crlevels,cslevels{1}) & strcmp(crnames,crname)) = 1;
%             end
            comp_Estimate(l) = rmat*re.Estimate + fmat*fe.Estimate;
            comp_pValue(l) = coefTest(LMEs{r}{t},fmat,0,'REContrast',rmat,'DFMethod','residual');
        end
        for cfe = 1:length(fe.Name)
            cfname = fe.Name{cfe};
            fprintf('%s %f %f\n', cfname, mean(comp_Estimate(strcmp(re.Group,"sub:chan") &...
                strcmp(re.Name,cfname))), fe.Estimate(cfe))
%             csubs = unique(re.Level(strcmp(re.Group,'sub')));
%             remeans =[];
%             for csubix = 1:length(csubs)
%                 csub = csubs{csubix};
%                 remeans(end+1) = mean(comp_Estimate(contains(re.Level,csub) & strcmp(re.Group,'sub:chan') & strcmp(re.Name,cfname)));                
%                 fprintf('%s %s %f %f\n', cfname, csub, remeans(end), comp_Estimate(strcmp(re.Level,csub) & strcmp(re.Group,'sub') & strcmp(re.Name,cfname)))
%             end
%             fprintf('%s %f %f\n', cfname, mean(remeans), fe.Estimate(cfe))
        end

        %%Compute additive estimates and p-values
%         comp_Estimate = nan(length(re),1);
%         comp_pValue = nan(length(re),1);
%         rlevels = unique(rnames.Level);
%         for f = 1:height(fnames)
%             cfname = fnames.Name{f};
%             fmat = zeros(1,height(fnames));
%             fmat(f) = 1;
%             for l = 1:length(rlevels)
%                 rmat = zeros(1,height(rnames));
%                 clevel = rlevels{l};
%                 cr = find(strcmp(rnames.Level,clevel) & strcmp(rnames.Name,cfname));
%                 cslevels = split(clevel,' ');
%                 rmat(1,cr) = 1;
%                 if length(cslevels) > 1
%                     rmat(strcmp(rnames.Level,cslevels{1}) & strcmp(rnames.Name,cfname)) = 1;
%                 end
%                 comp_Estimate(cr) = rmat*re + fmat*fe;
%                 comp_pValue(cr) = coefTest(LMEs{r}{t},fmat,0,'REContrast',rmat);%, 'DFMethod','residual');
%             end
%         end
        
        %%select channel random effect
        re = re(strcmp(re.Group,'sub:chan'),:);
        comp_Estimate = comp_Estimate(strcmp(re.Group,'sub:chan'));
        comp_pValue = comp_pValue(strcmp(re.Group,'sub:chan'));
        if t == 1
            chan_labels{r} = unique(re.Level);
            chan_coef{r} = NaN(size(model,2) + 1, length(chan_labels{r}), numel(LMEs{r}));
%             chan_lower{r} = chan_coef{r};
%             chan_upper{r} = chan_coef{r};
            chan_pval{r} = chan_coef{r};
        end
        
        chan_coef{r}(:,:,t) = reshape(comp_Estimate, size(chan_coef{r},1:2));
%         chan_lower{r}(:,:,t) = reshape(rndstats.Lower, size(chan_coef{r},1:2));
%         chan_upper{r}(:,:,t) = reshape(rndstats.Upper, size(chan_coef{r},1:2));
        chan_pval{r}(:,:,t) = reshape(comp_pValue, size(chan_coef{r},1:2));
    end
    chan_coef{r} = permute(chan_coef{r},[2,1,3]);
%     chan_lower{r} = permute(chan_lower{r},[2,1,3]);
%     chan_upper{r} = permute(chan_upper{r},[2,1,3]);
    chan_pval{r} = permute(chan_pval{r},[2,1,3]);
end

%% FDR correction
chan_qval = {};
for r = 1:numel(chan_pval)
    for ch = 1:size(chan_pval{r},1)
        for rg = 1:size(chan_pval{r},2)
            [~,~,~,chan_qval{r}(ch,rg,:)] = fdr_bh(chan_pval{r}(ch,rg,:));
        end
    end
end
%% Store in a structure
beta_chan =  [];
%conn_stats.LMEs = LMEs;
beta_chan.feature   = reg_lab;
beta_chan.time      = st.win_center;
beta_chan.coefs = chan_coef;
% beta_chan.lower = chan_lower;
% beta_chan.upper = chan_upper;
beta_chan.pvals = chan_pval;
beta_chan.qvals = chan_qval;
beta_chan.label = roi_list;
beta_chan.chan_label = chan_labels;
beta_chan.win_lim   = win_lim;
beta_chan.win_lim_s = hfa_time(win_lim);
beta_chan.dimord = 'chan_coef_time';

%% Estimate channel categories (Only if model EpnRPE_DifFV)
% if strcmp(model_id, 'EpnRPE_DifFB')
%     [cat_lab, ~, ~, ~, ~] = fn_puns_category_label_styles('puns');
%     for r = 1:numel(beta_chan.coefs)
%         
%         %array to store values
%         beta_chan.chancat_ix{r} = cell(4,1);
%         
%         % Find pRPE and nRPE channels with any significant time points
%         [pidx,~] = find(squeeze(beta_chan.qvals{r}(:,3,:) < .05));
%         [nidx,~] = find(squeeze(beta_chan.qvals{r}(:,4,:) < .05));
%         pidx = unique(pidx); nidx = unique(nidx);
%         
%         % find channels responding to either pRPE or nRPE (not both)
%         pchan = pidx(~ismember(pidx,nidx));
%         nchan = nidx(~ismember(nidx,pidx));
%         
%         % find channels responding to both pRPE and nRPE
%         common_chan = pidx(ismember(pidx,nidx));
%         
%         % calculate channels responding to uRPE (slnchan) and sRPE (rwdchan)
%         slnchan = [];
%         rwdchan = [];
%         for cch = 1:numel(common_chan)
%             % select current common channel
%             cchan = common_chan(cch);
%             
%             % get significant timepoints for pRPE and nRPE in this channel
%             psig = squeeze(beta_chan.qvals{r}(cchan,3,:) < .05);
%             nsig = squeeze(beta_chan.qvals{r}(cchan,4,:) < .05);
%             
%             % evaluate whether coefficients are +ve or -ve
%             ppos = squeeze(double(beta_chan.coefs{r}(cchan,3,:) > 0));
%             nneg = squeeze(double(beta_chan.coefs{r}(cchan,4,:) > 0));
%             
%             % if channel has only sig. positive coefs for pRPE and any sig.
%             % negative coef for nRPE (or vice versa), then classify as
%             % rwdchan (i.e. sRPE)
%             if (~ismember(0,ppos(psig)) & ismember(0,nneg(nsig))) |...
%                     (~ismember(0,nneg(nsig)) & ismember(0,ppos(psig)))
%                 rwdchan = [rwdchan;cchan];
%             else
%                 % if channel has any sig. negative coef for pRPE and any
%                 % sig. negative coef for nRPE, then clasify as sRPE ONLY
%                 % if not at same time points.
%                 
%                 %common significan times
%                 commtidx = find(double(psig).*double(nsig));
%                 
%                 %compare coefficients
%                 if sum(ppos(commtidx) ~= nneg(commtidx)) > 0
%                     rwdchan = [rwdchan;cchan];
%                 else
%                     % Everyhting else is uRPE
%                     slnchan = [slnchan;cchan];
%                 end
%             end
%         end
%         
%         cat_struc = [];
%         cat_struc.pRPE = pchan;
%         cat_struc.nRPE = nchan;
%         cat_struc.sRPE = rwdchan;
%         cat_struc.uRPE = slnchan;
%         
%         for ctg = 1:numel(cat_lab)
%             beta_chan.chancat_ix{r}{ctg} = cat_struc.(cat_lab{ctg});
%         end
%     end
%     beta_chan.chancat_label = cat_lab;
% end
%% Save
fprintf('=================== Saving channel effects ======================\n');
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
%out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef_REML.mat'];
save(out_fname,'-v7.3','beta_chan')
%save(out_fname,'beta_chan')
%% Old way of calculating significant channels
% chan_cat = cell(1,numel(beta_chan.coefs));
% for r = 1:numel(beta_chan.coefs)
%     % Find positive and negative channels with any significant time points
%        [allpidx,~] = find(squeeze(beta_chan.qvals{r}(:,3,:) < .05));
%        [allnidx,~] = find(squeeze(beta_chan.qvals{r}(:,4,:) < .05));
%         allpidx = unique(allpidx); allnidx = unique(allnidx);
%
%     % loop over time points
%     ntimes = size(beta_chan.coefs{r},3);
%     beta_chan.chancat_ix{r} = cell(4,ntimes);
%     for t = 1:ntimes
%         % Find positive and negative significant channels at current time
%         [pidx,~] = find(squeeze(beta_chan.qvals{r}(:,3,t) < .05));
%         [nidx,~] = find(squeeze(beta_chan.qvals{r}(:,4,t) < .05));
%         pidx = unique(pidx); nidx = unique(nidx);
%
%         % select channels that do not respond to opposite reward category
%         % in the whole trial
%         pchan = pidx(~ismember(pidx,allnidx));
%         nchan = nidx(~ismember(nidx,allpidx));
%
%         % get the sign of channel coefficient at current time
%         signcoef = double(beta_chan.coefs{r}(:, 3:4, t) >= 0);
%
%         % channels with equal coefficient sign (i.e. sRPE coding)
%         slnidx = find(signcoef(:,1) ~= signcoef(:,2));
%
%         % channels with different coefficient sign (i.e. uRPE coding)
%         rwdidx = find(signcoef(:,1) == signcoef(:,2));
%
%         % select significant uRPE channels
%         slnchan = pidx(ismember(pidx,nidx));
%         slnchan = slnchan(ismember(slnchan,slnidx));
%
%         % select non-significant
%         rwdchan = pidx(ismember(pidx,nidx));
%         rwdchan = rwdchan(ismember(rwdchan,rwdidx));
%
%         beta_chan.chancat_ix{r}{1,t} = pchan;
%         beta_chan.chancat_ix{r}{2,t} = nchan;
%         beta_chan.chancat_ix{r}{3,t} = slnchan;
%         beta_chan.chancat_ix{r}{4,t} = rwdchan;
%     end
% end
