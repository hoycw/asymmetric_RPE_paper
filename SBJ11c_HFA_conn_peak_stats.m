function SBJ11c_HFA_conn_peak_stats(proc_id, an_id, model_id, conn_id, stat_id, swap_Xcorr, rcn)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/conn/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir proc_id '_' model_id '_' an_id '_' stat_id '_' conn_id '_chancoef.mat'];
load(stats_fname)

if swap_Xcorr == 1
    for r = 1:numel(conn_stats_chan.coefs)
        conn_stats_chan.pair_label{r} = conn_stats_chan.pair_label{r}(1,[2,1]);
        for b = 1:size(conn_stats_chan.coefs{r},4)
            for ch = 1:size(conn_stats_chan.coefs{r},1)
                conn_stats_chan.coefs{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.coefs{r}(ch,:,:,b)));
                conn_stats_chan.pvals{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.pvals{r}(ch,:,:,b)));
                conn_stats_chan.qvals{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.qvals{r}(ch,:,:,b)));
%                 conn_stats_chan.lower{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.lower{r}(ch,:,:,b)));
%                 conn_stats_chan.upper{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.upper{r}(ch,:,:,b)));
                conn_stats_chan.peaks{r,b}.plat(ch) = conn_stats_chan.peaks{r,b}.plat(ch)*-1;
                conn_stats_chan.peaks{r,b}.nlat(ch) = conn_stats_chan.peaks{r,b}.nlat(ch)*-1;
                %conn_stats_chan.chan_label{r} = conn_stats_chan.chan_label{r}([2,1]);
            end
        end
    end
end
if ~isfield(conn_stats_chan,'nbins')
    conn_stats_chan.nbins = 1;
end
%% Get category data
d = conn_stats_chan.peaks{1,1};
cats = {'NS','pRPE','nRPE','sRPE','uRPE'};
d.subject = cellfun(@(x) x(1:4), d.pair, 'UniformOutput',false);
d.category_code = cell2mat(cellfun(@(x) find(strcmp(cats,x)), d.category, 'UniformOutput',false));
% Add categorical direction
d.ndir = sign(d.nlat);
d.pdir = sign(d.plat);
subs = unique(d.subject);

%% Get channel info
atlas_id  = 'Dx';
roi_id    = 'MPFCINS';%'gROI';%'mgROI';%'MPFCINS';%
rcn.hemi     = 'l';
rcn.mirror   = 1;
rcn.plot_roi = 'MPFCINS';
rcn = fn_process_recon_vars(rcn);
[elec_sbj, ~] = fn_load_grp_elec_ROI(subs,proc_id,atlas_id,roi_id,rcn);
elec = ft_appendsens([],elec_sbj{:});

%% Divide between anterior and posterior
for dix = 1:height(d)
    cMPFClab2 = split(d.pair{dix},'_to_');
    ch1 = cMPFClab2{1};
    SBJ = split(ch1,' ');   
    ch2 = [SBJ{1} ' ' cMPFClab2{2}];
    chix1 = find(strcmp(elec.label,ch1));
    chix2 = find(strcmp(elec.label,ch2));
    ycoord = elec.chanpos(chix1,2);
    if ycoord >= 12
        cant = 'ANT';
    else
        cant ='POST';
    end
    d.antpos{dix} = cant;
    d.x1{dix} = elec.chanpos(chix1,1);
    d.y1{dix} = elec.chanpos(chix1,2);
    d.z1{dix} = elec.chanpos(chix1,3);
    d.x2{dix} = elec.chanpos(chix2,1);
    d.y2{dix} = elec.chanpos(chix2,2);
    d.z2{dix} = elec.chanpos(chix2,3);
end
writetable(d, [stats_dir stat_id '_data_table.csv'], 'Delimiter',',')
%% Get category proportions
cat_prop = {};
for s = 1:numel(subs)
        props = zeros(5,1);
        cdata = d(strcmp(d.subject, subs{s}),:);
        numelec = max(4,length(unique(cdata.pair)));
        ctab = tabulate(cdata.category_code);
        if ~isempty(ctab)
            props(ctab(:,1)) = ctab(:,2) / numelec;
            cat_prop{end+1} = table(props);
        else
            cat_prop{end+1} = table(NaN(5,1));
        end
        cat_prop{end}.Properties.VariableNames = {'prop'};
        cat_prop{end}.subject = repmat({subs{s}}, length(props),1);
        cat_prop{end}.category_code = [0:4]';
end
cat_prop = vertcat(cat_prop{:});
writetable(cat_prop, [stats_dir stat_id '_pairs_categories_proportions.csv'])
%% Perform Kruskal Wallis test among categories
cat_prop2 = cat_prop(cat_prop.category_code > 0,:);
[~,ktbl,~] = kruskalwallis(cat_prop2.prop,cat_prop2.category_code, 'off');

%% report KW test
kstring = '\nKruskal-Wallis test among categories:\n\n';
kstring = [kstring sprintf('X2(%d) = %.2f, p = %.2f\n\n',ktbl{2,3}, ktbl{2,5}, ktbl{2,6})];
fprintf(kstring)
%% Perform Wilcoxon pairwise tests
wilcox_code = [4,1,2,3];
p=[]; h=[]; stats=[]; cpairs = []; cpair_labs = {};
for cat1 = 1:4
    for cat2 =1:4 
        if cat2 > cat1
            code1 = wilcox_code(cat1);
            code2 = wilcox_code(cat2);
            x = cat_prop2.prop(cat_prop2.category_code == code1);
            y = cat_prop2.prop(cat_prop2.category_code == code2);
            [p(end+1), h(end+1), cstat] = signrank(x,y);
            stats(end+1) = cstat.signedrank;
            cpairs(end+1,:) = [code1,code2];
            cpair_labs{end+1} = [cats{code1+1} '-' cats{code2+1}];
        end
    end
end
cat_stats = [];
cat_stats.p = p;
[~,~,~,cat_stats.q] = fdr_bh(p);
cat_stats.h = h;
cat_stats.stat = stats;
cat_stats.pairs = cpairs;
cat_stats.pair_labels = cpair_labs;
%% Report
cstring = '\nPairwise Wilcoxon signed-rank test between categories:\n\n';
for cat = 1:numel(cat_stats.pair_labels)
    curpair = cat_stats.pair_labels{cat};
    cstring = [cstring  sprintf('%s: statistic = %02.2f, p = %.3f, q = %.3f\n',...
                curpair, cat_stats.stat(1,cat),cat_stats.p(cat),cat_stats.q(cat))];
end
cstring = [cstring '\n'];
fprintf(cstring)

%% Calculate medians and IQR

dstring = '\nMedians and inter-quartile ranges (percent of pairs):\n\n';
for cat = 1:4
    cprops = cat_prop.prop(cat_prop.category_code == cat);
    mdn = median(cprops,'omitnan')*100;
    ciqr = prctile(cprops,[25,75])*100;
    dstring = [dstring sprintf('%s: MDN = %.2f, IQR = %.2f - %.2f\n',cats{cat+1},mdn,ciqr(1),ciqr(2))];
end
fprintf(dstring)

%% Proportion of inhibitory channels
coefs = {'p','n','p','p'};
istring = '\n';
for c = 2:numel(cats)
    ccat = cats{c};
    n1 = sum(strcmp(d.category,ccat) & sign(d.([coefs{c-1} 'peak']))==-1); 
    n2 = sum(strcmp(d.category,ccat));
    istring = [istring sprintf('Proportion of inhibitory %s peaks: %d / %d = %.2f\n',ccat,n1,n2,n1/n2)];
end
fprintf(istring)
%% Lag stats
%% Number of positive vs negative lags
cat_excs = {[0,2],[0,1]};
ctypes = {'p','n'};
cnames = {'positive','negative'};
formulas = {'lat ~ 1 + (1|subject)';
            'lat ~ category  + (1+category|subject)';
            'lat ~ category  + antpos + (1+category + antpos|subject)';
            'lat ~ category*antpos + (1+category*antpos|subject)'};
anat_models = [];
anat_string = '\nTest of peak lags according to category and anatomy\n';
for ct = 1:length(ctypes)
    anat_string = [anat_string, '\n' cnames{ct},' coefficients\n'];
    ctype = ctypes{ct};
    ccatix = ~ismember(d.category_code,cat_excs{ct});
    for f = 1:length(formulas)
        cmodel = ['m' num2str(f)];
        anat_models.(ctype).(cmodel) = fitlme(d(ccatix,:),[ctype,formulas{f}]);
        if f > 1
            cLRT = compare(anat_models.(ctype).(['m' num2str(f-1)]),...
                           anat_models.(ctype).(cmodel));
            anat_string = sprintf('%s\nnull: %s\nalternative: %s\n\n',anat_string,formulas{f-1},formulas{f});
            anat_string = sprintf('%sdeltaDF: %d\nLRstat: %f\np: %f\n',anat_string,cLRT.deltaDF(2),cLRT.LRStat(2),cLRT.pValue(2));
        end
    end
end
fprintf(anat_string)
%% Perform Kruskal Wallis test among categories
pcatix = ~ismember(d.category_code,[0,2]);
ncatix = ~ismember(d.category_code,[0,1]);
[~,ptbl,~] = kruskalwallis(d.plat(pcatix), d.category_code(pcatix), 'off');
[~,ntbl,~] = kruskalwallis(d.nlat(ncatix), d.category_code(ncatix), 'off');
%% report KW test
pstring = '\nKruskal-Wallis test among categories (positive RPE coefficient):\n\n';
pstring = [pstring sprintf('X2(%d) = %.2f, p = %.2f\n\n',ptbl{2,3}, ptbl{2,5}, ptbl{2,6})];
nstring = '\nKruskal-Wallis test among categories (negative RPE coefficient):\n\n';
nstring = [nstring sprintf('X2(%d) = %.2f, p = %.2f\n\n',ntbl{2,3}, ntbl{2,5}, ntbl{2,6})];
fprintf(pstring);fprintf(nstring);
%% Perform Wilcoxon pairwise tests for negative coefficients
wilcox_code = [4,1,2,3];
p=[]; h=[]; stats=[]; cpairs = []; cpair_labs = {};
for cat1 = 2:4
    for cat2 =2:4 
        if cat2 > cat1
            code1 = wilcox_code(cat1);
            code2 = wilcox_code(cat2);
            x = d.nlat(d.category_code == code1);
            y = d.nlat(d.category_code == code2);
            [p(end+1), h(end+1), cstat] = ranksum(x,y);
            stats(end+1) = cstat.ranksum;
            cpairs(end+1,:) = [code1,code2];
            cpair_labs{end+1} = [cats{code1+2} '-' cats{code2+2}];
        end
    end
end
lag_stats = [];
lag_stats.p = p;
[~,~,~,lag_stats.q] = fdr_bh(p);
lag_stats.h = h;
lag_stats.stat = stats;
lag_stats.pairs = cpairs;
lag_stats.pair_labels = cpair_labs;
%% Report
lstring = '\nPairwise Wilcoxon rank-sum test between categories for peak lags (n coeff):\n\n';
for cat = 1:numel(lag_stats.pair_labels)
    curpair = lag_stats.pair_labels{cat};
    lstring = [lstring  sprintf('%s: statistic = %02.2f, p = %.3f, q = %.3f\n',...
                curpair, lag_stats.stat(1,cat),lag_stats.p(cat),lag_stats.q(cat))];
end
lstring = [lstring '\n'];
fprintf(lstring)
%% Perform Wilcoxon pairwise tests for positive coefficients
wilcox_code = [4,1,2,3];
p=[]; h=[]; stats=[]; cpairs = []; cpair_labs = {};
for cat1 = 2:4
    for cat2 =2:4
        if cat2 > cat1
            code1 = wilcox_code(cat1);
            code2 = wilcox_code(cat2);
            x = d.plat(d.category_code == code1);
            y = d.plat(d.category_code == code2);
            [p(end+1), h(end+1), cstat] = ranksum(x,y);
            stats(end+1) = cstat.ranksum;
            cpairs(end+1,:) = [code1,code2];
            cpair_labs{end+1} = [cats{code1+2} '-' cats{code2+2}];
        end
    end
end
lag_stats = [];
lag_stats.p = p;
[~,~,~,lag_stats.q] = fdr_bh(p);
lag_stats.h = h;
lag_stats.stat = stats;
lag_stats.pairs = cpairs;
lag_stats.pair_labels = cpair_labs;
%% Report
lpstring = '\nPairwise Wilcoxon rank-sum test between categories for peak lags (p coeff):\n\n';
for cat = 1:numel(lag_stats.pair_labels)
    curpair = lag_stats.pair_labels{cat};
    lpstring =[lpstring  sprintf('%s: statistic = %02.2f, p = %.3f, q = %.3f\n',...
                curpair, lag_stats.stat(1,cat),lag_stats.p(cat),lag_stats.q(cat))];
end
lpstring = [lpstring '\n'];
fprintf(lpstring)

%% Calculate medians and IQR of peak lats
coefs = {'p','n'};
cats2 = {{'pRPE','sRPE','uRPE'},{'nRPE','sRPE','uRPE'}};
mstring = '';
for c = 1:numel(coefs)
    mstring = [mstring sprintf('\nMedians and inter-quartile ranges of %s coefficient peak latencies:\n\n',coefs{c})];
    ccoef = [coefs{c} 'lat'];
    ccats = cats2{c};
    for cat = 1:numel(ccats)
        clats = d.(ccoef)(strcmp(d.category,ccats{cat}));
        mdn = median(clats,'omitnan');
        ciqr = prctile(clats,[25,75]);
        mstring = [mstring sprintf('%s: MDN = %.2f, IQR = %.2f - %.2f\n',ccats{cat},mdn,ciqr(1),ciqr(2))];
    end
end
fprintf(mstring)

%% Save reports
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
stats_file = sprintf('%s%s_%s_%s_conn_categories_stats.txt', stats_dir,...
                    model_id, an_id, stat_id);
disp('saving report')
fid = fopen(stats_file,'w');
fprintf(fid,[kstring, cstring, dstring,istring,anat_string,...
             pstring,nstring,lstring,lpstring, mstring]);
fclose(fid);
