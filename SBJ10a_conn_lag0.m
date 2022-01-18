function SBJ10a_conn_lag0(SBJ, method, proc_id, an_id, roi_id, atlas_id, Nbins, chansel)
%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else; root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Data
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);

% Load elec, HFA data and active electrodes
load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'.mat']);
load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'_actv_mn100.mat']);

% Check if more than one frequency, error for now
if numel(hfa.freq)>1
    error('HFA has more than one frequency, can''t run on that for now!');
end

%%
if strcmp(method,'corr')
    Nbins = 1;
end
%% Select channels in ROIs
[roi_list,~, roi_field] = fn_roi_label_styles(roi_id);

if strcmp(chansel, 'maxactive')
    r_ix = [];
    %select max active channel per ROI:
    for r = 1:numel(roi_list)
        cr_ix = find(strcmp(elec.(roi_field), roi_list{r}));
        [minqval,midx] = min(min(actv.qvals(cr_ix,:),[],2)); %select channel with lowest pval
        if minqval < .05
            r_ix = [r_ix;cr_ix(midx)]; % make sure it is under .05
        end
    end
elseif ~strcmp(chansel, 'maxactive') | numel(r_ix) < numel(roi_list)
    r_ix = [];
    %else select all channels
    for r = 1:numel(roi_list)
        r_ix = [r_ix; find(strcmp(elec.(roi_field),roi_list{r}))];
    end
    if strcmp(chansel,'active')
        %select active channels if required
        r_ix = r_ix(ismember(r_ix, find(actv.actv_ch == 1)));
    end
end
roi_labs = elec.(roi_field)(r_ix);
chan_labs = elec.label(r_ix);
%%  Compute pairwise connectivity between channels in ROIs:
fs = (length(hfa.time) - 1)/ (hfa.time(end) - hfa.time(1));

conn = []; connidx = 0;
for rl1 = 1:length(roi_list)
    for rl2 = 1:length(roi_list)
        if rl2 > rl1
            connidx = connidx + 1;
            chidx1 = find(strcmp(roi_labs, roi_list{rl1}));
            chidx2 = find(strcmp(roi_labs, roi_list{rl2}));
            conn.pair_labels{connidx}{1} = roi_list{rl1};
            conn.pair_labels{connidx}{2} = roi_list{rl2};
            conn.chann_labels{connidx}{1} = chan_labs(chidx1);
            conn.chann_labels{connidx}{2} = chan_labs(chidx2);
            conn.coeff{1,connidx} = NaN(size(hfa.powspctrm,1),numel(Nbins),numel(chidx1),numel(chidx2));
            for tr = 1:size(hfa.powspctrm,1)
                fprintf('calculating 0-lag connectivity for subject %s trial %d / %d\n',SBJ,tr, size(hfa.powspctrm,1))
                for n = 1:numel(Nbins)
                    for ch1 = 1:numel(chidx1)
                        for ch2 = 1:numel(chidx2)
                            c1 = squeeze(hfa.powspctrm(tr,chidx1(ch1),:,:));
                            c2 = squeeze(hfa.powspctrm(tr,chidx2(ch2),:,:));
                            if strcmp(method,'MI')
                            c1bin = rh_discretize(c1, Nbins(n), 'UniCount');
                            c2bin = rh_discretize(c2, Nbins(n), 'UniCount');
                            conn.coeff{1,connidx}(tr,n,ch1,ch2) = MutualInformation(c1bin, c2bin);
                            clear c1bin c2bin
                            elseif strcmp(method,'corr')
                               curcoefs =  corrcoef(c1, c2);
                               conn.coeff{1,connidx}(tr,n,ch1,ch2) = curcoefs(1,2);
                            else
                             error('Method not known (specify either "MI" or "corr")')        
                            end
                            clear c1 c2
                        end
                    end
                end
            end
        end
    end
end
%% Complete structure and save
conn.label = chan_labs;
conn.roi = roi_labs;
conn.dimord = 'rpt_bin_chan_chan';
conn.fsampleorig = fs;
conn.trialinfo = hfa.trialinfo;
conn.nbins = Nbins;
out_fname = [SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'_',method,'_0lag','.mat'];
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',out_fname);
fprintf('===================================================\n');
save(out_fname,'conn','-v7.3');
%% Make figures reporting individual connectivity
fig_dir = [root_dir 'PRJ_Error/results/conn/' SBJ '/' proc_id '_' an_id '/' ];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

for pc = 1:length(conn.coeff)
    fig = figure;
    set(gcf, 'PaperType','a2')
    set(gcf, 'PaperOrientation','landscape')
    for n = 1:numel(Nbins)
        cpdata = squeeze(mean(conn.coeff{1,pc}(:,n,:,:),1));
        subplot(ceil(length(Nbins)/3),3,n)
        imagesc(cpdata); colorbar();
        title(num2str(Nbins(n)))
        set(gca,'XTick',1:1:size(cpdata,2),'XTickLabel',conn.chann_labels{pc}{2})
        set(gca,'YTick',1:1:size(cpdata,1),'YTickLabel',flipud(conn.chann_labels{pc}{1}))
    end
    sgtitle(sprintf('0-lag %s - %s %s',...
        conn.pair_labels{pc}{1},conn.pair_labels{pc}{2},method))
    
    figname = sprintf('%s%s_%s_0-lagged_%s_%s',...
        fig_dir, SBJ,method,conn.pair_labels{pc}{1},conn.pair_labels{pc}{2});
    print(fig, figname, '-dpdf','-fillpage')
    close(fig)
end
end
