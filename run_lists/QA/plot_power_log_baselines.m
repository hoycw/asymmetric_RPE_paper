% function plot_power_log_baselines(SBJ,proc_id,an_id)
% Filter, extract power, plot distribution of power values

%% Settings
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
SBJ       = 'CP24';
proc_id   = 'main_ft';
ch_ix     = 1;
bsln_type = 'zscore';
save_fig  = 0;
fig_ftype = 'png';

fig_dir = [root_dir 'PRJ_Error/results/log_power_testing/' SBJ '/'];% an_id '/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

%% Set up paths
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Compute parameters
an.std_thresh = 3;
an.fois = [2:2:12 14:4:30 40:10:200];%2.^[2:1/4:8]; % Center frequencies
an.octave = 0.25;
an.bplim = zeros([numel(an.fois) 2]);
for f = 1:numel(an.fois)
    an.bplim(f,1) = 2^(-an.octave/2)*an.fois(f);
    an.bplim(f,2) = 2^(an.octave/2)*an.fois(f);
end
% an.bplim = round(an.bplim);

% an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
an.evnt_lab    = 'S';           % event around which to cut trials
% trial_lim_s will NOT be full of data! the first and last t_ftimwin/2 epochs will be NaNs
an.trial_lim_s = [-0.25 3.01];      % window in SEC for cutting trials
an.demean_yn   = 'no';             % z-score for HFA instead
an.bsln_evnt   = 'S';
% an.bsln_type   = 'zboot';
an.bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction
an.bsln_boots  = 500;              % repetitions for non-parametric stats

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));

% Toss sampleinfo which can mess up trial cutting
if isfield(data,'sampleinfo')
    data = rmfield(data,'sampleinfo');
end

%% Select Channel(s)
cfgs = [];
cfgs.channel = data.label(1);%an.
roi = ft_selectdata(cfgs,data);
% roi = data;

%% Compute TFRs
cfg = [];
cfg.hilbert  = 'abs';
cfg.bpfilter = 'yes';
cfg.bpfreq   = [];      % to be filled by looping through foi_center
cfg.bpfiltord = 2;
cfg.channel  = 'all';

pow = cell(size(an.fois));
maxs = zeros(size(an.fois));
for f = 1:numel(an.fois)
    cfg.bpfreq = an.bplim(f,:);%[an.fois(f_ix)-bws(f_ix)/2 an.fois(f_ix)+bws(f_ix)/2];
    fprintf('\n------> filtering: %.03f - %.03f\n', cfg.bpfreq(1), cfg.bpfreq(2));
    pow{f} = ft_preprocessing(cfg,roi);
    maxs(f) = max(pow{f}.trial{1});
end

%% Cut into trials
bsln_events = trl_info.trl_onset;
if strcmp(an.evnt_lab,'S')
    % Check that baseline will be included in data cut to trial_lim_s
    if an.trial_lim_s(1) < an.bsln_lim(1)
        error(['ERROR: an.trial_lim_s does not include an.bsln_lim for an_id = ' an_id]);
    end
    % Cut to desired trial_lim_s
    cut_lim = round(an.trial_lim_s*roi.fsample);
elseif strcmp(an.evnt_lab,'R')
    % Check that baseline will be included in data cut to trial_lim_s
    if an.trial_lim_s(1)+min(trl_info.response_time) < an.bsln_lim(1)
        error(['ERROR: an.trial_lim_s does not include an.bsln_lim for an_id = ' an_id]);
    end
    % Cut to max_RT+trial_lim_s(2) to include S baseline + full R-locked trial_lim_s
    max_RT  = max(trl_info.rt);
    cut_lim = round([an.trial_lim_s(1) max_RT+an.trial_lim_s(2)]*roi.fsample);
elseif strcmp(an.evnt_lab,'F')
    % Cut from S baseline to full F-locked:
    %   [trl_onset+an.trial_lim_s(1) to fb+an.trial_lim_s(2)]
    cut_lim = round([an.trial_lim_s(1) trl_info.prdm.trl_len+an.trial_lim_s(2)]*roi.fsample);
else
    error(['Unknown an.evnt_lab: ' an.evnt_lab]);
end

pow_trl  = cell(size(pow));
lpow_trl = cell(size(pow));
for f = 1:numel(an.fois)
    pow_trl{f} = fn_ft_cut_trials_equal_len(pow{f},bsln_events,ones(size(trl_info.trl_n)),cut_lim);
    % Log Transform
    tmp = pow{f};
    tmp.trial{1} = log(tmp.trial{1});
    lpow_trl{f} = fn_ft_cut_trials_equal_len(tmp,bsln_events,ones(size(trl_info.trl_n)),cut_lim);
end

%% Baseline Correction
bsln_pow  = cell(size(pow));
bsln_lpow = cell(size(pow));
for f = 1:numel(an.fois)
    switch bsln_type
        case {'zscore', 'zboot', 'demean', 'my_relchange'}
            bsln_pow{f} = fn_bsln_ft_filtered(pow_trl{f},an.bsln_lim,bsln_type,an.bsln_boots);
            bsln_lpow{f} = fn_bsln_ft_filtered(lpow_trl{f},an.bsln_lim,bsln_type,an.bsln_boots);
        case {'db','relchange'}
            cfgbsln = [];
            cfgbsln.baseline     = an.bsln_lim;
            cfgbsln.baselinetype = bsln_type;
            cfgbsln.parameter    = 'trial';
            bsln_pow{f}  = ft_freqbaseline(cfgbsln,pow_trl{f});
            bsln_lpow{f} = ft_freqbaseline(cfgbsln,lpow_trl{f});
        otherwise
            error(['No baseline implemented for an.bsln_type: ' bsln_type]);
    end
end

%% Grab trial matrices
bsln_pow_mat  = cell(size(pow));
bsln_lpow_mat = cell(size(pow));
for f = 1:numel(an.fois)
    bsln_pow_mat{f}  = nan([numel(trl_info.trl_n) size(bsln_pow{f}.trial{1},2)]);
    bsln_lpow_mat{f} = nan([numel(trl_info.trl_n) size(bsln_lpow{f}.trial{1},2)]);
    for t = 1:numel(trl_info.trl_n)
        bsln_pow_mat{f}(t,:)  = bsln_pow{f}.trial{t};
        bsln_lpow_mat{f}(t,:) = bsln_lpow{f}.trial{t};
    end
end

%% Plot Raw Distributions
bins = linspace(0,max(maxs),500);
fig_name = [SBJ '_bsln_power_' roi.label{1}];
figure('Name',fig_name,'Units','Normalized','outerposition',[0 0 1 1]);
n_rows = ceil((numel(an.fois)+1)/2);
subplot(n_rows,2,1);

% Plot original data
histogram(roi.trial{1});
title(['orig: skew = ' num2str(skewness(roi.trial{1})) '; kurt = ' num2str(kurtosis(roi.trial{1}))]);

means = zeros(size(an.fois));
vars  = zeros(size(an.fois));
skews = zeros(size(an.fois));
kurts = zeros(size(an.fois));
n_out = zeros(size(an.fois));
for f = 1:numel(an.fois)
    means(f) = mean(bsln_pow_mat{f}(:));
    vars(f)  = std(bsln_pow_mat{f}(:));
    skews(f) = skewness(bsln_pow_mat{f}(:));
    kurts(f) = kurtosis(bsln_pow_mat{f}(:));
    thresh = means(f) + vars(f)*an.std_thresh;
    n_out(f) = 100*sum(bsln_pow_mat{f}(:)>thresh)/numel(bsln_pow_mat{f}(:));
    if f<=n_rows-1; c=1; else; c=2;end
    sp = fn_rowcol2subplot_ix(n_rows,2,mod(f,n_rows)+1,c);
    subplot(n_rows,2,sp);
    histogram(bsln_pow_mat{f}(:),bins);
    line([thresh thresh],ylim,'Color','r');
    title([num2str(an.bplim(f,1),3) '-' num2str(an.bplim(f,2),3) ' Hz: skew = ' num2str(skews(f),'%.2f')...
        '; kurt = ' num2str(kurts(f),'%.2f') '; %out = ' num2str(n_out(f),'%.2f')]);
end
% suptitle('Raw Power Distributions');

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Log Transforms
bins = linspace(-10,5,500);
fig_name = [SBJ '_log_bsln_power_' roi.label{1}];
figure('Name',fig_name,'Units','Normalized','outerposition',[0 0 1 1]);
subplot(n_rows,2,1);

% Plot original data
histogram(roi.trial{1});
title(['original: skew = ' num2str(skewness(roi.trial{1})) '; kurt = ' num2str(kurtosis(roi.trial{1}))]);

lmeans = zeros(size(an.fois));
lvars  = zeros(size(an.fois));
lskews = zeros(size(an.fois));
lkurts = zeros(size(an.fois));
ln_out = zeros(size(an.fois));
for f = 1:numel(an.fois)
    lmeans(f) = mean(bsln_lpow_mat{f}(:));
    lvars(f)  = std(bsln_lpow_mat{f}(:));
    lskews(f) = skewness(bsln_lpow_mat{f}(:));
    lkurts(f) = kurtosis(bsln_lpow_mat{f}(:));
    thresh = lmeans(f) - lvars(f)*an.std_thresh;
    ln_out(f) = 100*sum(bsln_lpow_mat{f}(:)<thresh)/numel(pow{f}.trial{1});
    if f<=n_rows-1; c=1; else; c=2;end
    sp = fn_rowcol2subplot_ix(n_rows,2,mod(f,n_rows)+1,c);
    subplot(n_rows,2,sp);
    histogram(bsln_lpow_mat{f}(:),bins);
    line([thresh thresh],ylim,'Color','r');
    title([num2str(an.bplim(f,1),3) '-' num2str(an.bplim(f,2),3) ' Hz: skew = ' num2str(skews(f),'%.2f')...
        '; kurt = ' num2str(kurts(f),'%.2f') '; %out = ' num2str(ln_out(f),'%.2f')]);
end

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Statistical Moments across Frequencies
fig_name = [SBJ '_RawLog_power_summary_' roi.label{1}];
figure('Name',fig_name,'Units','Normalized','outerposition',[0 0 0.6 0.8]);

% Power Spectrum
subplot(3,2,1, 'XScale', 'log', 'YScale', 'log'); hold on;
[fft_data,freqs] = pwelch(roi.trial{1},2048,0,2048,roi.fsample);
plot(freqs,fft_data);
xlim([1 350]);
ax = gca;
ax.XTick = [4 8 12 25 30 60 80 100 120 180 200 240 300 360 420];
xlabel('frequency');
ylabel('Power');
title('PSD');

% Outliers
subplot(3,2,2); hold on;
plot(an.fois,n_out,'r');
plot(an.fois,ln_out,'b');
xlabel('frequency');
ylabel('% Outliers');
legend('raw','log');
title(['Gaussian Outliers (>' num2str(std_thresh) ' SDs)']);

% Means
subplot(3,2,3);hold on;
plot(an.fois,means,'r');
plot(an.fois,lmeans,'b');
xlabel('frequency');
ylabel('Mean');
legend('raw','log');
scatter(0,mean(roi.trial{1}),25,'k');
title('Mean');

% Standard Deviations
subplot(3,2,4);hold on;
plot(an.fois,vars,'r');
plot(an.fois,lvars,'b');
xlabel('frequency');
ylabel('StD');
legend('raw','log');
%scatter(0,std(roi.trial{1}),25,'k');
title('Standard Deviation');

% Skewness
subplot(3,2,5);hold on;
plot(an.fois,skews,'r');
plot(an.fois,lskews,'b');
xlabel('frequency');
ylabel('Skew');
legend('raw','log');
scatter(0,skewness(roi.trial{1}),25,'k');
title('Skewness');

% Kurtosis
subplot(3,2,6);hold on;
plot(an.fois,kurts,'r');
plot(an.fois,lkurts,'b');
xlabel('frequency');
ylabel('Kurtosis');
legend('raw','log');
scatter(0,kurtosis(roi.trial{1}),25,'k');
title('Kurtosis');

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end


% end