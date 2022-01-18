function SBJ11a_HFA_conn_mGLM(SBJ, proc_id, an_id, model_id, stat_id, conn_id, pmask, nbins)
%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Data
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
% eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
if ~strcmp(st.an_style,'mGLM'); error('This script is only for mass GLM!'); end

% Load behavioral data, model, and connectivity data
load([SBJ_vars.dirs.events,SBJ,'_bhv_',proc_id,'_final.mat']);
load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'_', conn_id,'.mat']);


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

% trimm data
if isfield(conn,'time')
    tridx = find(conn.time >= st.stat_lim(1) & conn.time <= st.stat_lim(2));
    conn.time = conn.time(tridx);
    for cmp = 1:length(conn.coeff)
        conn.coeff{cmp} = conn.coeff{cmp}(:,:,:,:,tridx);
    end
else
    conn.time = 0;
    st.stat_lim = [0,0];
end

%% Check and select bins
if isfield(conn,'nbins')
   bidx = find(ismember(conn.nbins, nbins));
   if ~isempty(bidx)
       disp(['slecting nbins ' num2str(conn.nbins(bidx))]);
       conn.nbins= conn.nbins(bidx);
       for cf = 1:numel(conn.coeff)
           conn.coeff{cf} = conn.coeff{cf}(:,bidx,:,:,:); 
       end
   else
       disp(['could not find nbins ' num2str(nbins)])
       disp(['selecting all bins available: ' num2str(conn.nbins)])
   end
else
    conn.nbins = 1;
end
%% Run GLM
% reshape data for GLM
if ~any(strcmp(reg_lab,'offset'))
    reg_lab = ['offset' reg_lab];
    model   = [ones(size(model,1), 1) model];
end

% Regress rt ifasked for
if st.regress_rt
    if ~any(strcmp(reg_lab,'rt'))
        disp('RT will be regressed out')
        model  = [model, zscore(bhv.rt)];
        reg_lab = [reg_lab 'rt'];
    end
end

beta = conn;
beta.model     = model;
beta.feature   = reg_lab;
beta.win_lim   = st.stat_lim;
beta.trial = {};
beta.boot = {};
beta.zscore = {};
beta.bootmean = {};
beta.bootstd = {};
beta.pval = {};
beta.qval = {};
beta.dimord = strrep(beta.dimord, 'rpt', 'feat');
beta = rmfield(beta,'coeff');
beta = rmfield(beta,'trialinfo');

for cmp = 1:length(conn.coeff)
    
    Y = reshape(conn.coeff{cmp},size(conn.coeff{cmp},1),[],size(conn.coeff{cmp},5));
    
    % Bootstrap results
    
    %beta.max_hfa_z = max_z;
    % w2.trial     = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
    %beta.boot      = zeros([numel(beta.feature) length(hfa.label) length(beta.time) st.n_boots]);
    % w2.pval      = w2.trial;
    %beta.win_lim_s = hfa.time(win_lim);
    
    % Compute GLM for real model
    beta.trial{cmp} = fn_mass_GLM(model,Y,0);
    beta.trial{cmp}= reshape(beta.trial{cmp},[size(beta.trial{cmp},1),size(conn.coeff{cmp},2:5)]);
    beta.boot{cmp} = zeros([size(beta.trial{cmp},1:5),st.n_boots]);
    
    % Permuting labels %
    fprintf('boot #: ');
    for boot_ix = 1:st.n_boots
        %     m = sprintf(' permutation %d/%d', boot_ix, n_boots);
        %     fprintf([b m]); b = repmat('\b',[1 length(m)]);
        fprintf('%i..',boot_ix);
        rand_model = model(randperm(size(model,1)),:);
        beta.boot{cmp}(:,:,:,:,:,boot_ix) = reshape(fn_mass_GLM(rand_model,Y,0),size(beta.trial{cmp},1:5));
        if mod(boot_ix,20)==0
            fprintf('\n');
        end
    end
    
    % compute pvals and do multiple comparison corrections
    beta.pval{cmp} = sum(bsxfun(@ge,abs(beta.boot{cmp}-median(beta.boot{cmp},6)),...
        abs(beta.trial{cmp}-median(beta.boot{cmp},6))),6)./st.n_boots;
    
    % Multiple Comparisons Correction within Channel pair
    beta.qval{cmp} = nan(size(beta.pval{cmp},1:5));
    for ch_ix1 = 1:numel(beta.chann_labels{cmp}{1})
        for ch_ix2 = 1:numel(beta.chann_labels{cmp}{2})
            [~, ~, ~, beta.qval{cmp}(:,:,ch_ix1,ch_ix2,:)] = fdr_bh(squeeze(beta.pval{cmp}(:,:,ch_ix1,ch_ix2,:)));
        end
    end
    
    % Z-score beta statistics for later group use
    beta.zscore{cmp}   = norminv(1-beta.pval{cmp},0,1);
    beta.bootmean{cmp} = mean(beta.boot{cmp},6);
    beta.bootstd{cmp}  = std(beta.boot{cmp},[],6);
    beta.zscore{cmp}(isinf(beta.zscore{cmp})) = norminv(1-1/st.n_boots/2,0,1);
end
%% Plot results
fig_dir = [root_dir 'PRJ_Error/results/conn/' SBJ '/' proc_id '_' an_id '/' ];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

xticks =  1 : floor(length(beta.time) / 4) : length(beta.time);
xticklabels = beta.time(xticks);
maxncols = 4;
for cc = 1:numel(beta.pair_labels)
    roi1 = beta.pair_labels{cc}{1};
    roi2 = beta.pair_labels{cc}{2};
    cur_chans1 = beta.chann_labels{cc}{1};
    cur_chans2 = beta.chann_labels{cc}{2};
    npcols = min(numel(cur_chans1), maxncols);
    if npcols == 1; npcols = 2; end
    for ff = 1:numel(beta.feature)
        for n = 1:numel(beta.nbins)
            fig = figure;
            set(gcf, 'PaperType','a2')
            set(gcf, 'PaperOrientation','landscape')
            for cch = 1:numel(cur_chans1)
                subplot(ceil(numel(cur_chans1)/npcols),npcols,cch)
                cpdata = squeeze(beta.trial{cc}(ff,n,cch,:,:));
                if pmask == 1
                    cpdata = cpdata.*(squeeze(beta.qval{cc}(ff,n,cch,:,:))<=0.05);
                endfig_dir = [root_dir 'PRJ_Error/results/HFA/'
                imagesc(cpdata); colorbar();
                title(cur_chans1{cch})
                set(gca, 'XTick',xticks,'XTickLabel',xticklabels)
                if cch == 1
                    set(gca,'YTick',1:size(cur_chans2),'YTickLabel',flipud(cur_chans2))
                else
                    set(gca,'YTickLabel',[])
                end
                sgtitle(sprintf('connectivity from %s to %s - %s', roi1, roi2, beta.feature{ff}))
            end
            figname = [fig_dir, SBJ, '_mGLM_conn_', beta.pair_labels{cc}{1},'_to_',...
                      beta.pair_labels{cc}{2},'_',beta.feature{ff},'_' conn_id,...
                      '_' num2str(beta.nbins(n)),'bins'];
            print(fig, figname, '-dpdf','-fillpage')
            close(fig)
        end
    end
end

% fig_dir = [root_dir 'PRJ_Error/results/conn/' SBJ '/' proc_id '_' an_id '/' ];
% for cmp = 1:numel(beta.trial)
%     for n = 1:numel(beta.nbins)
%         fig = figure;
%         set(gcf, 'PaperType','a2')
%         if numel(beta.time) > 1
%             set(gcf, 'PaperOrientation','portrait')
%         else
%             set(gcf, 'PaperOrientation','landscape')
%         end
%         for f = 1:numel(beta.feature)
%             cur_beta = squeeze(beta.trial{cmp}(f,n,:,:,:));
%             cur_qval = squeeze(beta.qval{cmp}(f,n,:,:,:));
%             if numel(beta.time) > 1
%                 cur_beta = reshape(cur_beta,[],size(cur_beta,3));
%                 cur_qval = reshape(cur_qval,[],size(cur_qval,3));
%                 xlabs = beta.time;
%                 ylabs = {};
%                 for ch1 = 1:numel(beta.chann_labels{cmp}{1})
%                     for ch2 = 1:numel(beta.chann_labels{cmp}{2})
%                         ylabs(end+1) = {[beta.chann_labels{cmp}{1}{ch1},' to ',...
%                                          beta.chann_labels{cmp}{2}{ch2}]};
%                     end
%                 end
%                 yticks = 1:numel(ylabs);
%                 xticks = [1,ceil(numel(beta.time)/2),numel(beta.time)];
%                 xlabs = beta.time(xticks);
%             else
%                 ylabs = beta.chann_labels{cmp}{1};
%                 yticks = 1:numel(ylabs);
%                 xlabs = beta.chann_labels{cmp}{2};
%                 xticks = 1:numel(xlabs);
%             end
%             subplot(1,numel(beta.feature),f)
%             imagesc(cur_beta.*(cur_qval<.05));colorbar()
%             set(gca, 'XTick', xticks, 'XTickLabel', xlabs)
%             set(gca, 'YTick', yticks,'YTickLabel', flipud(ylabs))
%             set(gca, 'FontSize', 6)
%             title(beta.feature{f})
%         end
%         sgtitle(sprintf('connectivity from %s to %s', beta.pair_labels{cmp}{1},...
%                         beta.pair_labels{cmp}{2}))
%         figname = [fig_dir, SBJ, '_mGLM_conn_', beta.pair_labels{cmp}{1},'_to_',...
%                    beta.pair_labels{cmp}{2},'_', conn_id, '_' num2str(beta.nbins(n)),'bins'];
%         print(fig, figname, '-dpdf','-fillpage')
%         close(fig)
%     end
% end
%% Save Results
out_fname = [SBJ_vars.dirs.stats SBJ '_mGLM_conn_' model_id '_' stat_id '_' an_id '_' conn_id '.mat'];
save(out_fname,'-v7.3','beta');
