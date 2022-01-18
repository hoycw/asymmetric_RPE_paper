function SBJ10c_conn_win(SBJ, method, proc_id, an_id, roi_id, atlas_id, twin, tlags, stepsize, Nbins, chansel)
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
%%  Compute pairwise MI connectivity between channels in ROIs:
fs = (length(hfa.time) - 1)/ (hfa.time(end) - hfa.time(1));
newfs = 1 / stepsize;
tvec = tlags(1) : stepsize : tlags(2);

conn = []; connidx = 0;
for rl1 = 1:length(roi_list)
    for rl2 = 1:length(roi_list)
        if rl1 ~= rl2
            connidx = connidx + 1;
            chidx1 = find(strcmp(roi_labs, roi_list{rl1}));
            chidx2 = find(strcmp(roi_labs, roi_list{rl2}));
            conn.pair_labels{connidx}{1} = roi_list{rl1};
            conn.pair_labels{connidx}{2} = roi_list{rl2};
            conn.chann_labels{connidx}{1} = chan_labs(chidx1);
            conn.chann_labels{connidx}{2} = chan_labs(chidx2);
            conn.coeff{1,connidx} = NaN(size(hfa.powspctrm,1),numel(Nbins),numel(chidx1),numel(chidx2),numel(tvec));
            for tr = 1:size(hfa.powspctrm,1)
                fprintf('calculating connectivity for subject %s trial %d / %d - %s to %s\n',...
                    SBJ,tr, size(hfa.powspctrm,1),roi_list{rl1},roi_list{rl2})
                for n = 1:numel(Nbins)
                    for ch1 = 1:numel(chidx1)
                        for ch2 = 1:numel(chidx2)
                            for t = 1:numel(tvec)
                                tidx0 = nearest(hfa.time,twin(1)) : nearest(hfa.time,twin(2));
                                tidx1 = nearest(hfa.time,twin(1)+tvec(t)) : nearest(hfa.time,twin(2)+tvec(t));
                                c1 = squeeze(hfa.powspctrm(tr,chidx1(ch1),:,tidx0));
                                c2 = squeeze(hfa.powspctrm(tr,chidx2(ch2),:,tidx1));
                                if size(c1) == size(c2)
                                    if strcmp(method, 'MI')
                                        c1bin = rh_discretize(c1, Nbins(n), 'UniCount');
                                        c2bin = rh_discretize(c2, Nbins(n), 'UniCount');
                                        conn.coeff{1, connidx}(tr,n,ch1,ch2,t) = MutualInformation(c1bin, c2bin);
                                        clear c1bin c2bin
                                    elseif strcmp(method, 'corr')
                                        curcoefs =  corrcoef(c1, c2);
                                        conn.coeff{1,connidx}(tr,n,ch1,ch2,t) = curcoefs(1,2);
                                    else
                                        error('Method not known (specify either "MI" or "corr")')
                                    end
                                    clear c1 c2
                                else
                                    error('Time window mismatch. Please check your time window and lags')
                                end
                            end
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
conn.dimord = 'rpt_bin_chan_chan_time';
conn.fsampleorig = fs;
conn.time = tvec;
conn.trialinfo = hfa.trialinfo;
conn.nbins = Nbins;
out_fname = [SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'_',method,'_win.mat'];
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',out_fname);
fprintf('===================================================\n');
save(out_fname,'conn','-v7.3');

%% Make figures reporting individual connectivity
fig_dir = [root_dir 'PRJ_Error/results/conn/' SBJ '/' proc_id '_' an_id '/' ];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

xticks =  1 : floor(length(tvec) / 4) : length(tvec);
xticklabels = tvec(xticks);
maxncols = 4;
for cc = 1:numel(conn.pair_labels)
    roi1 = conn.pair_labels{cc}{1};
    roi2 = conn.pair_labels{cc}{2};
    cur_chans1 = conn.chann_labels{cc}{1};
    cur_chans2 = conn.chann_labels{cc}{2};
    npcols = min(numel(cur_chans1), maxncols);
    if npcols == 1; npcols = 2; end
    for n = 1:numel(Nbins)
        fig = figure;
        set(gcf, 'PaperType','a2')
        set(gcf, 'PaperOrientation','landscape')
        for cch = 1:numel(cur_chans1)
            subplot(ceil(numel(cur_chans1)/npcols),npcols,cch)
            cpdata = squeeze(mean(conn.coeff{cc}(:,n,cch,:,:),1));
            if size(squeeze(cpdata),2) == 1
                plot(tvec,cpdata)
                set(gca, 'XTick',xticklabels,'XTickLabel',xticklabels)
                title([cur_chans1{cch} ' to ' cur_chans2{1}])
            else
                imagesc(cpdata); colorbar();
                title(cur_chans1{cch})
                set(gca, 'XTick',xticks,'XTickLabel',xticklabels)
                if cch == 1
                    set(gca,'YTick',1:size(cur_chans2),'YTickLabel',flipud(cur_chans2))
                else
                    set(gca,'YTickLabel',[])
                end
            end
            sgtitle(sprintf('connectivity from %s to %s', roi1, roi2))
        end
        figname = sprintf('%s%s_%s_twin%d_to_%d_%s_to_%s_%dbins_win',...
            fig_dir, SBJ, method, twin(1)*1000, twin(2)*1000, roi1 ,roi2, Nbins(n));
        print(fig, figname, '-dpdf','-fillpage')
        close(fig)
    end
end
