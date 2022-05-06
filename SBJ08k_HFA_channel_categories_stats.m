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
tables = {};
for r = 1:numel(beta_chan.chancat_ix)   
    tables{r} = table(beta_chan.chan_label{r});
    tables{r}.Properties.VariableNames = {'channel'};
    tables{r}.subject = cellfun(@(x) x(1:4), tables{r}.channel, 'UniformOutput',false);
    tables{r}.region = repmat({beta_chan.label{r}},height(tables{r}),1);
    tables{r}.region_code = repmat([r-1],height(tables{r}),1);
    tables{r}.category = repmat({'nonsignificant'}, height(tables{r}),1);
    tables{r}.category_code = repmat([0], height(tables{r}),1);
    for cat = 1:numel(beta_chan.chancat_label)
        tables{r}.category(beta_chan.chancat_ix{r}{cat}) = repmat(beta_chan.chancat_label(cat),...
                                                                  numel(beta_chan.chancat_ix{r}{cat}),1);
        tables{r}.category_code(beta_chan.chancat_ix{r}{cat}) = repmat(cat,...
                                                                  numel(beta_chan.chancat_ix{r}{cat}),1);                                                             
    end
end
d = vertcat(tables{:});
writetable(d, [stats_dir 'channel_categories_' beta_chan.pval_type,'vals.csv'])
%% stats
d2 = d(d.category_code > 0,:);
[~, ~, multinom_stats] = mnrfit(double(d2.region_code), double(d2.category_code));

%% Calculate proportion data
subs = unique(d.subject);
rois = unique(d.region);

cat_prop = {};
for s = 2:numel(subs)
    for r = 1:numel(rois)
        props = zeros(5,1);
        cdata = d(strcmp(d.subject, subs{s}) & strcmp(d.region, rois{r}),:);
        numelec = max(4,length(unique(cdata.channel)));
        ctab = tabulate(cdata.category_code);
        props(ctab(:,1)+1) = ctab(:,2) / numelec;
        cat_prop{end+1} = table(props);
        cat_prop{end}.Properties.VariableNames = {'prop'};
        cat_prop{end}.subject = repmat({subs{s}}, length(props),1);
        cat_prop{end}.region = repmat({rois{r}}, length(props),1);
        cat_prop{end}.category_code = [0:4]';
        size(cat_prop{end})
    end
end
cat_prop = vertcat(cat_prop{:});

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

%% Perform Wilcoxon test between categories
p=[]; h=[]; stats=[]; cpairs = []; cpair_labs = {};
for cat1 = 1:4
    for cat2 = 1:4
        if cat2 > cat1 
            x = cat_prop.prop(cat_prop.category_code == cat1);
            y = cat_prop.prop(cat_prop.category_code == cat2);
            [p(end+1), h(end+1), cstat] = signrank(x,y);
            stats(end+1) = cstat.signedrank;
            cpairs(end+1,:) = [cat1,cat2];
            cpair_labs{end+1} = [cat_lab{cat1} '-' cat_lab{cat2}];
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