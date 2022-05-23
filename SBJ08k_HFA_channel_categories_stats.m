function SBJ08k_HFA_channel_categories_stats(an_id, model_id, stat_id)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname)

%% estimate categories
[cat_lab, ~, ~, ~, ~] = fn_puns_category_label_styles('puns');

%% Put together in a table
mnrfit_code = [2,3,4,1];
%mnrfit_code = [1,4,3,2];
tables = {};
for r = 1:numel(beta_chan.chancat_ix)   
    tables{r} = table(beta_chan.chan_label{r});
    tables{r}.Properties.VariableNames = {'channel'};
    tables{r}.subject = cellfun(@(x) x(1:4), tables{r}.channel, 'UniformOutput',false);
    tables{r}.region = repmat({beta_chan.label{r}},height(tables{r}),1);
    tables{r}.region_code = repmat([r-1],height(tables{r}),1);
    tables{r}.category = repmat({'nonsignificant'}, height(tables{r}),1);
    tables{r}.category_code = repmat([0], height(tables{r}),1);
    tables{r}.mnrfit_code = tables{r}.category_code;
    for cat = 1:numel(beta_chan.chancat_label)
        nchan_cat = numel(beta_chan.chancat_ix{r}{cat});
        tables{r}.category(beta_chan.chancat_ix{r}{cat}) = repmat(beta_chan.chancat_label(cat),nchan_cat,1);
        tables{r}.category_code(beta_chan.chancat_ix{r}{cat}) = repmat(cat,nchan_cat,1);
        tables{r}.mnrfit_code(beta_chan.chancat_ix{r}{cat}) = repmat(mnrfit_code(cat),nchan_cat,1);
    end
end
d = vertcat(tables{:});
writetable(d, [stats_dir 'channel_categories_' beta_chan.pval_type,'vals.csv'])

%% Calculate proportion data
subs = unique(d.subject);
rois = unique(d.region);

cat_prop = {};
for s = 1:numel(subs)
    for r = 1:numel(rois)
        props = zeros(5,1);
        cdata = d(strcmp(d.subject, subs{s}) & strcmp(d.region, rois{r}),:);
        numelec = max(4,length(unique(cdata.channel)));
        ctab = tabulate(cdata.category_code);
        if ~isempty(ctab)
            props(ctab(:,1)+1) = ctab(:,2) / numelec;
            cat_prop{end+1} = table(props);
        else
            cat_prop{end+1} = table(NaN(5,1));
        end
        cat_prop{end}.Properties.VariableNames = {'prop'};
        cat_prop{end}.subject = repmat({subs{s}}, length(props),1);
        cat_prop{end}.region = repmat({rois{r}}, length(props),1);
        cat_prop{end}.category_code = [0:4]';
    end
end
cat_prop = vertcat(cat_prop{:});
writetable(cat_prop, [stats_dir 'channel_categories_proportions_' beta_chan.pval_type,'vals.csv'])
%% stats
d2 = d(d.mnrfit_code > 0,:);
[~, ~, mnstats] = mnrfit(double(d2.region_code), double(d2.mnrfit_code));
%% report
mstring = '\nMultinomial test of channel categories proportions:\n\nIntercept\n\n';
refcat = cat_lab{mnrfit_code==4}; 
for cfs = 1:size(mnstats.beta, 2)
    curcat = cat_lab{mnrfit_code==cfs};
    mstring = [mstring  sprintf('%s - %s: odds ratio = %.2f, CI(95) = [%.2f %.2f], p = %.2f\n',...
                curcat, refcat, exp(mnstats.beta(1,cfs)),...
                exp(mnstats.beta(1,cfs)-1.96*(mnstats.se(1,cfs))),...
                exp(mnstats.beta(1,cfs)+1.96*(mnstats.se(1,cfs))),...
                mnstats.p(1,cfs))];
end
mstring = [mstring '\nRegion\n\n'];
for cfs = 1:size(mnstats.beta, 2)
    curcat = cat_lab{mnrfit_code==cfs};
    mstring = [mstring  sprintf('%s - %s: odds ratio = %.2f, CI(95) = [%.2f %.2f], p = %.2f\n',...
                curcat, refcat, exp(mnstats.beta(2,cfs)),...
                exp(mnstats.beta(1,cfs)-1.96*(mnstats.se(2,cfs))),...
                exp(mnstats.beta(1,cfs)+1.96*(mnstats.se(2,cfs))),...
                mnstats.p(2,cfs))];
end
mstring = [mstring '\n'];
fprintf(mstring)
%% Perform Wilcoxon test between regions
p=NaN(1,4); h=NaN(1,4); stats=NaN(1,4);
for cat = 1:4
    x = cat_prop.prop(cat_prop.category_code == cat & strcmp(cat_prop.region, rois{1}));
    y = cat_prop.prop(cat_prop.category_code == cat & strcmp(cat_prop.region, rois{2}));
    %[p(cat), h(cat), stats(cat)] = signrank(x,y);
    [p(cat), h(cat), cstat] = signrank(x,y);
    stats(cat) = cstat.signedrank;
end
roi_stats = [];
roi_stats.p = p;
[~,~,~,roi_stats.q] = fdr_bh(p);
roi_stats.h = h;
roi_stats.stat = stats;
roi_stats.category_codes = 1:4;
roi_stats.categories = cat_lab;
roi_stats.rois = rois;
%% Report
rstring = sprintf('\nPairwise Wilcoxon signed-rank test of categories between regions (%s - %s):\n\n',rois{1},rois{2});
for cat = 1:numel(roi_stats.categories)
    curcat = roi_stats.categories{cat};
    rstring = [rstring  sprintf('%s: statistic = %02.2f, p = %.2f, q = %.2f\n',...
                curcat, roi_stats.stat(1,cat),roi_stats.p(cat),roi_stats.q(cat))];
end
rstring = [rstring '\n'];
fprintf(rstring)
%% Perform Kruskal Wallis test among categories
cat_prop2 = cat_prop(strcmp(cat_prop.region, 'INS'),:);
cat_prop2.prop = mean([cat_prop.prop(strcmp(cat_prop.region, 'INS')),...
                       cat_prop.prop(strcmp(cat_prop.region, 'MPFC'))],2);
cat_prop2 = cat_prop2(cat_prop2.category_code > 0,:);
[~,ktbl,~] = kruskalwallis(cat_prop2.prop,cat_prop2.category_code, 'off');

%% report KW test
kstring = '\nKruskal-Wallis test among categories across regions:\n\n';
kstring = [kstring sprintf('X2(%d) = %.2f, p = %.2f\n\n',ktbl{2,3}, ktbl{2,5}, ktbl{2,6})];
fprintf(kstring)

%% Perform Wilcoxon pairwise tests
wilcox_code = [3,4,1,2];
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
            cpair_labs{end+1} = [cat_lab{code1} '-' cat_lab{code2}];
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
cstring = '\nPairwise Wilcoxon signed-rank test between categories across regions:\n\n';
for cat = 1:numel(cat_stats.pair_labels)
    curpair = cat_stats.pair_labels{cat};
    cstring = [cstring  sprintf('%s: statistic = %02.2f, p = %.3f, q = %.3f\n',...
                curpair, cat_stats.stat(1,cat),cat_stats.p(cat),cat_stats.q(cat))];
end
cstring = [cstring '\n'];
fprintf(cstring)
%% Calculate medians and IQR
dstring = '\nMedians and inter-quartile ranges (percent of channels):\n';

for r = 1:numel(rois)
    roi = rois{r};
    dstring = [dstring sprintf('\nRegion: %s\n\n', roi)];
    for cat = 1:4
        cprops = cat_prop.prop(cat_prop.category_code == cat & strcmp(cat_prop.region,roi));
        mdn = median(cprops,'omitnan')*100;
        ciqr = prctile(cprops,[25,75])*100;
        dstring = [dstring sprintf('%s: MDN = %.2f, IQR = %.2f - %.2f\n',cat_lab{cat},mdn,ciqr(1),ciqr(2))];
    end
end
fprintf(dstring)

%% Calculate number and percentage of inhibitory channels
reg_ix = [3,4,3,3];
istring = '\n\nProportion of channels with "inhibitory" responsivenes:\n';
for cat = 1:numel(beta_chan.chancat_label)
    creg = reg_ix(cat);
    rProps = []; rCount = []; rN = [];
    clabel = beta_chan.chancat_label{cat};
    istring = [istring '\n' clabel '\n\n'];
    for r = 1:numel(beta_chan.chancat_ix)
        cat_sign = [];
        rlabel = beta_chan.label{r};
        cchans = beta_chan.chancat_ix{r}{cat};
        for ch = 1:numel(cchans)
            cch = cchans(ch);
            sig = beta_chan.([beta_chan.pval_type 'vals']){r}(cch,creg,:) < .05;
            cbetas = squeeze(beta_chan.coefs{r}(cch,creg,:) .* sig);
            pks = [unique(min(cbetas)), unique(max(cbetas))];
            [~, pix] = max(abs(pks));
            cat_sign(end + 1) = sign(pks(pix));           
        end
        rN(end+1) = length(cat_sign);
        rCount(end+1) = sum(cat_sign == -1);
        rProps(end+1) = rCount(end) / rN(end);
        istring = [istring sprintf('%s: %d / %d = %.2f\n',rlabel,...
                                   rCount(end),rN(end),rProps(end))];
    end
    istring = [istring, sprintf('Total: %d / %d = %.2f\n', sum(rCount),...
                                sum(rN),sum(rCount)/sum(rN))];
end
fprintf(istring)
%% Save reports
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
stats_file = sprintf('%s%s_%s_%s_hfa_categories_stats.txt', stats_dir,...
                    model_id, an_id, stat_id);
fid = fopen(stats_file,'w');
fprintf(fid,[mstring,rstring,kstring,cstring, dstring, istring]);
fclose(fid);