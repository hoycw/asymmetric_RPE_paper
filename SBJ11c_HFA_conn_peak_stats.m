function SBJ11c_HFA_conn_peak_stats(proc_id, an_id, model_id, conn_id, swap_Xcorr)

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
stats_fname = [stats_dir proc_id '_' model_id '_' an_id '_' conn_id '_chancoef.mat'];
load(stats_fname)

if swap_Xcorr == 1
    for r = 1:numel(conn_stats_chan.coefs)
        conn_stats_chan.pair_label{r} = conn_stats_chan.pair_label{r}(1,[2,1]);
        for b = 1:size(conn_stats_chan.coefs{r},4)
            for ch = 1:size(conn_stats_chan.coefs{r},1)
                conn_stats_chan.coefs{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.coefs{r}(ch,:,:,b)));
                conn_stats_chan.pvals{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.pvals{r}(ch,:,:,b)));
                conn_stats_chan.qvals{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.qvals{r}(ch,:,:,b)));
                conn_stats_chan.lower{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.lower{r}(ch,:,:,b)));
                conn_stats_chan.upper{r}(ch,:,:,b) = fliplr(squeeze(conn_stats_chan.upper{r}(ch,:,:,b)));
                %conn_stats_chan.chan_label{r} = conn_stats_chan.chan_label{r}([2,1]);
            end
        end
    end
end
if ~isfield(conn_stats_chan,'nbins')
    conn_stats_chan.nbins = 1;
end
%% find peaks
% ntimes = numel(conn_stats_chan.time);
peaks = cell(numel(conn_stats_chan.coefs),size(conn_stats_chan.coefs{1},4));
lags = peaks;
sig_ix = peaks;
for r = 1:numel(conn_stats_chan.coefs)
    for b = 1:size(conn_stats_chan.coefs{r},4)
        % findpeaks
        peaks{r,b} = NaN(size(conn_stats_chan.coefs{r},1),2);
        lags{r,b} = peaks{r,b};
        for ch =  1:length(peaks{r,b})
            for creg = 3:4
                cchdata = squeeze(conn_stats_chan.coefs{r}(ch,creg,:,b));
                [cpeak, latix] = findpeaks(abs(cchdata));
                [~,pidx] = max(cpeak);
                peaks{r, b}(ch,creg-2) = cchdata(latix(pidx(1)));
                lags{r, b}(ch,creg-2) =  conn_stats_chan.time(latix(pidx(1)));
            end
        end
    end
    sig_ix{r,b} = mean(conn_stats_chan.qvals{r}(:,3:4,:,b) < 0.05, 3) > 0;
end

%% Do stats and report
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
stats_file = sprintf('%s%s_%s_%s_%s_hfa_peak_lags_stats_%s_to_%s.txt', stats_dir,...
                proc_id, model_id, an_id, conn_id, conn_stats_chan.pair_label{r}{1},...
                conn_stats_chan.pair_label{r}{2});
fid = fopen(stats_file,'w');
for r = 1:size(lags,1)
    for b = 1:size(lags,2)
        nbins = conn_stats_chan.nbins(b);
        stitle = sprintf('Tests of peak connectivity lags from %s to %s (nbis = %d):\n\n',...
                        conn_stats_chan.pair_label{r}{1},conn_stats_chan.pair_label{r}{2},nbins);
        [~,p,ci,stats] = ttest(lags{r,b}(:,1), lags{r,b}(:,2));
        [~,p2,ci2,stats2] = ttest2(lags{r,b}(:,1), lags{r,b}(:,2));
        [~,p3,ci3,stats3] = ttest2(lags{r,b}(sig_ix{r,b}(:,1),1),  lags{r,b}(sig_ix{r,b}(:,2),2));
        
        % report
        line1 = ['paired t-test between pRPE and nRPE peak lags (all channel pairs): ',...
            sprintf('\n\nt(%d) = %f\nM = %f s\nCI = [%f %f]\np = %f\n\n',...
            stats.df,stats.tstat, mean(diff(fliplr(lags{r,b})')), ci(1), ci(2), p)];
        line2 = ['two-samples t-test between pRPE and nRPE peak lags (all channel pairs): ',...
            sprintf('\n\nt(%d) = %f\nM = %f s\nCI = [%f %f]\np = %f\n\n',...
            stats2.df,stats2.tstat, diff(fliplr(mean(lags{r,b},1))), ci2(1), ci2(2), p2)];
        line3 = ['two-samples t-test between pRPE and nRPE peak lags (significant channel pairs): ',...
            sprintf('\n\nt(%d) = %f\nM = %f s\nCI = [%f %f]\np = %f\n\n',...
            stats3.df, stats3.tstat,...
            mean(lags{r,b}(sig_ix{r,b}(:,1),1)) - mean(lags{r,b}(sig_ix{r,b}(:,2),2)),...
            ci3(1), ci3(2), p3)];
        
        disp([stitle, line1, line2, line3]);
        fprintf(fid,[stitle, line1, line2, line3])
    end
end
fclose(fid);
%% Save to text file
