function SBJ11b_HFA_conn_mLME(SBJs, proc_id, an_id, model_id, conn_id, stat_id, lme_formula,lags)
%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults
fig_dir = [root_dir 'PRJ_Error/results/conn/GRP/'];
%% Load Data for each subject and stack into tables:
conn_data = [];
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('loading subject %s\n',SBJ)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
    % eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
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
                            conn_data{pl}{s,cnt,cb}{dcount,3} = [lab1 '_to_' lab2];
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
            if st.rank_transform == 1
                disp('transforming to rankit')
                ctable.y = norminv((tiedrank(ctable.y)-.5) / height(ctable))*std(ctable.y)+mean(ctable.y);
            end
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
resid = {};
fittedv = {};
% lmresid={];e_formula = ['y~sRPE + uRPE + EV + (1 + sRPE + uRPE + EV | sub) +', ...
%                 '(1 + sRPE + uRPE + EV | sub:chan)'];
clix = find(time>=lags(1) & time<=lags(2));
for pl = 1:numel(conn_tables)
    nb = size(conn_tables{pl},2);
    nt = length(clix); %size(conn_tables{pl},1);
    LMEs{pl} = cell(nt,nb);
    pvals{pl} = NaN(size(model,2) + 1,nt,nb);
    coefs{pl} = NaN(size(model,2) + 1,nt,nb);
    lower{pl} = NaN(size(model,2) + 1,nt,nb);
    upper{pl} = NaN(size(model,2) + 1,nt,nb);
    resid{pl} = NaN(height(conn_tables{pl}{1,1}),nt,nb);
    fittedv{pl} = NaN(height(conn_tables{pl}{1,1}),nt,nb);
    for b = 1:nb
        for t = 1:nt
            ct = clix(t);
            %             fprintf('normalizing data\n')
            %             conn_tables{pl}{t,b}.y = (conn_tables{pl}{t,b}.y - mean(conn_tables{pl}{t,b}.y))...
            %                                         /  std(conn_tables{pl}{t,b}.y);
            fprintf('fitting label pair %d time point %d and bin %d\n',pl,t,b)
            LMEs{pl}{t,b} = fitlme(conn_tables{pl}{ct,b}, lme_formula);
            pvals{pl}(:,t,b) = LMEs{pl}{ct,b}.Coefficients.pValue;
            coefs{pl}(:,t,b) = LMEs{pl}{ct,b}.Coefficients.Estimate;
            lower{pl}(:,t,b) = LMEs{pl}{ct,b}.Coefficients.Lower;
            upper{pl}(:,t,b) = LMEs{pl}{ct,b}.Coefficients.Upper;
            resid{pl}(:,t,b) = residuals(LMEs{pl}{ct,b});
            fittedv{pl}(:,t,b) = fitted(LMEs{pl}{ct,b});
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
conn_stats.time = time(clix);
conn_stats.dimord = 'coef_time_bin';

%% Plot residuals qqplot
for pl = 1:numel(conn_tables)
    for b = 1:nb
        cresid = squeeze(resid{pl}(:,:,b));
        [~,lpval,~,~] = lillietest(reshape(cresid,[numel(cresid),1]));
        
        cf = figure('units','normalized','outerposition',[0 0 1 1],...
            'PaperOrientation','Landscape');
        subplot(1,2,1)
        qqplot(reshape(cresid,[numel(cresid),1])/std(cresid,[],'all'))
        title(sprintf('Lillefors pval = %.3f',lpval))
        %ylim([-5,5])

        subplot(1,2,2)
        [f, x] = hist(cresid, 50); % Create histogram from a normal distribution.
        g = 1 / sqrt(2 * pi) * exp(-0.5 * (x/std(cresid,[],'all') ) .^ 2); % pdf of the normal distribution

        % METHOD 1: DIVIDE BY SUM
        bar(x, f / trapz(x, f)); hold on
        plot(x, g / trapz(x, g), 'r'); hold off  
        xlim([-5*std(cresid,[],'all'), 5*std(cresid,[],'all')])
        %histfit(reshape(cresid,[numel(cresid),1]))
        title(sprintf('median = %.3f',median(cresid,'all')))
        
        plot_fname1 = [fig_dir proc_id '_' model_id '_' an_id '_' stat_id...
                      '_' conn_id '_residuals_qqplot_' pair_labels{pl}{1}...
                      '_to_' pair_labels{pl}{2} '_' num2str(b) '.pdf'];
        %
        print(plot_fname1,cf,'-dpdf','-fillpage')
    end
end
%close all
%% Residuals vs fitted
for pl = 1:numel(conn_tables)
    for b = 1:nb
        cf = figure('units','normalized','outerposition',[0 0 1 1],...
            'PaperOrientation','Landscape');
        cresid = resid{pl}(:,:,b);
        cfitted = fittedv{pl}(:,:,b);
        s = scatter(reshape(cfitted,[numel(cfitted),1]),...
            reshape(cresid,[numel(cresid),1]),'.k');
        s.MarkerEdgeAlpha = .05;
        s.MarkerFaceAlpha = .05;
        title([pair_labels{pl}{1} ' to ' pair_labels{pl}{2}])
        xlabel('fitted')
        ylabel('residuals')
        plot_fname2 = [fig_dir proc_id '_' model_id '_' an_id '_' stat_id...
                      '_' conn_id '_residuals_vs_fitted_' pair_labels{pl}{1}...
                      '_to_' pair_labels{pl}{2} '_' num2str(b) '.png'];
        print(plot_fname2,cf,'-dpng')%,'-fillpage')
    end
end
%%
close all
%% Save
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir proc_id '_' model_id '_' an_id '_' stat_id '_' conn_id '.mat'];
save(out_fname,'-v7.3','conn_stats')
%% Extract electrode-wise coefficients.
chan_coef = {};
% chan_lower = {};
% chan_upper = {};
chan_pval = {};
chan_labels = cell(size(LMEs));
for r = 1:numel(LMEs)
%     chan_labels{r} = unique(rnames.Level);
%     chan_coef{r} = NaN(size(model,2) + 1, length(chan_labels{r}),...
%                          size(LMEs{r},1),size(LMEs{r},2));
% %                 chan_lower{r} = chan_coef{r};
% %                 chan_upper{r} = chan_coef{r};
%     chan_pval{r} = chan_coef{r};
    for t = 1:size(LMEs{r},1)
        for b = 1:size(LMEs{r},2)
            %[~,~,rndstats] =  randomEffects(LMEs{r}{t,b});
            
            [re,rnames,~] =  randomEffects(LMEs{r}{t,b},'DFMethod','residual');%'satterthwaite')%;
            [fe,fnames,~] = fixedEffects(LMEs{r}{t,b});
            fprintf('estimating random effects for pair %d, timepoint %d and bin %d\n',r,t,b)
            nr = length(re);
            nf = length(fe);
            
            %%Compute additive estimates and p-values
            comp_Estimate = nan(nr,1);
            comp_pValue = nan(nr,1);
            %rlevels = unique(rnames.Level);
%             for f = 1:height(fnames)
%                 cfname = fnames.Name{f};
%                 disp(cfname)
%            fmat = zeros(1,height(fnames));
%                 fmat(f) = 1;
            cfnames=fnames.Name;
            crnames=rnames.Name;
            crlevels=rnames.Level;
            for l = 1:nr
                crname = crnames{l};
                fmat = zeros(1,nf);
                rmat = zeros(1,nr);
                fmat(1,strcmp(cfnames,crname)) = 1;
                rmat(1,l) = 1;
                clevel = crlevels{l};
                %disp(clevel)
%                cslevels = split(clevel,' ');
%                 if length(cslevels) > 1
%                     rmat(1,strcmp(crlevels,cslevels{1}) & strcmp(crnames,crname)) = 1;
%                 end
                comp_Estimate(l) = rmat*re + fmat*fe;
                comp_pValue(l) = coefTest(LMEs{r}{t},fmat,0,'REContrast',rmat,'DFMethod','residual');
            end
            for cfe = 1:length(fe)
                cfname = fnames.Name{cfe};
                fprintf('%s %f %f\n', cfname, mean(comp_Estimate(strcmp(rnames.Group,"sub:chan") &...
                    strcmp(rnames.Name,cfname))), fe(cfe))
            end
            %end

            %%select channel random effect
            rnames = rnames(strcmp(rnames.Group,'sub:chan'),:);
            comp_Estimate = comp_Estimate(strcmp(rnames.Group,'sub:chan'));
            comp_pValue = comp_pValue(strcmp(rnames.Group,'sub:chan'));
        
            %select channel random effect
            %rndstats = rndstats(strcmp(rndstats.Group,'sub:chan'),:);
            if t == 1 & b == 1
                cchan_labels = unique(rnames.Level);
                cchan_coef = NaN(size(model,2) + 1, length(cchan_labels),...
                    size(LMEs{r},1),size(LMEs{r},2));
    %                 chan_lower{r} = chan_coef{r};
    %                 chan_upper{r} = chan_coef{r};
                cchan_pval = cchan_coef;
            end
            cchan_coef(:,:,t,b) = reshape(comp_Estimate,size(cchan_coef,1:2));
%             chan_lower{r}(:,:,t,b) = reshape(rndstats.Lower,size(chan_coef{r},1:2));
%             chan_upper{r}(:,:,t,b) = reshape(rndstats.Upper,size(chan_coef{r},1:2));
            cchan_pval(:,:,t,b) = reshape(comp_pValue,size(cchan_coef,1:2));
        end
    end
    chan_coef{r} = permute(cchan_coef,[2,1,3,4]);
    chan_labels{r} = cchan_labels;
%     chan_lower{r} = permute(chan_lower{r},[2,1,3,4]);
%     chan_upper{r} = permute(chan_upper{r},[2,1,3,4]);
    chan_pval{r} = permute(cchan_pval,[2,1,3,4]);
end

%% FDR correction
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

%% Store in a structure
conn_stats_chan =  [];
%conn_stats.LMEs = LMEs;
conn_stats_chan.feature   = reg_lab;
conn_stats_chan.time      = time(clix);
conn_stats_chan.coefs = chan_coef;
% conn_stats_chan.lower = chan_lower;
% conn_stats_chan.upper = chan_upper;
conn_stats_chan.pvals = chan_pval;
conn_stats_chan.qvals = chan_qval;
conn_stats_chan.pair_label = pair_labels;
conn_stats_chan.chan_label = chan_labels;
conn_stats_chan.dimord = 'chan_coef_time_bin';

%% % Estimate channel categories
%[cat_lab, ~, ~, ~, ~] = fn_puns_category_label_styles('puns');
%
% for r = 1:numel(conn_stats_chan.coefs)
%
%     %array to store values
%     conn_stats_chan.chancat_ix{r} = cell(4,size(conn_stats.coefs{r},4));
%     for b = 1:size(conn_stats_chan.coefs{r},4)
%         % Find pRPE and nRPE channels with any significant time points
%         [pidx,~] = find(squeeze(conn_stats_chan.qvals{r}(:,3,:,b) < .05));
%         [nidx,~] = find(squeeze(conn_stats_chan.qvals{r}(:,4,:,b) < .05));
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
%             psig = squeeze(conn_stats_chan.qvals{r}(cchan,3,:,b) < .05);
%             nsig = squeeze(conn_stats_chan.qvals{r}(cchan,4,:,b) < .05);
%
%             % evaluate whether coefficients are +ve or -ve
%             ppos = squeeze(double(conn_stats_chan.coefs{r}(cchan,3,:,b) > 0));
%             nneg = squeeze(double(conn_stats_chan.coefs{r}(cchan,4,:,b) > 0));
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
%                 %common significant times
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
%             conn_stats_chan.chancat_ix{r}{ctg,b} = cat_struc.(cat_lab{ctg});
%         end
%     end
% end
%
% conn_stats_chan.chancat_label = cat_lab;

%% save
fprintf('=================== Saving channel effects ======================\n');
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir proc_id '_' model_id '_' an_id '_' stat_id '_' conn_id '_chancoef.mat'];
save(out_fname,'-v7.3','conn_stats_chan')
