function plot_power_RawVsLog(SBJ,proc_id,an_id)
% Filter, extract power, plot distribution of power values

% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));

% Toss sampleinfo which can mess up trial cutting
if isfield(data,'sampleinfo')
    data = rmfield(data,'sampleinfo');
end

%% Select Channel(s)
cfgs = [];
cfgs.channel = data.label(1);
roi = ft_selectdata(cfgs,data);
%clear data;

%% Plot PSD
figure;
[fft_data,freqs] = pwelch(roi.trial{1},2048,0,2048,roi.fsample);
loglog(freqs,fft_data);
xlim([1 350]);
ax = gca;
ax.XTick = [4 8 12 25 30 60 80 100 120 180 200 240 300 360 420];

%% Compute TFRs
% all cfg_tfr options are specified in the an_vars
fois = [8:2:20 25:5:60 70:10:200];%2.^[2:1/4:8]; % Center frequencies
%bws  = [repmat(2,[1 10]) repmat(5,[1 8]) repmat(10,[1 14])];
octave = 0.25;
bplim = zeros([numel(fois) 2]);
for f = 1:numel(fois)
    bplim(f,1) = 2^(-octave/2)*fois(f);
    bplim(f,2) = 2^(octave/2)*fois(f);
end

%%
cfg_hfa = [];
cfg_hfa.hilbert  = 'abs';
cfg_hfa.bpfilter = 'yes';
cfg_hfa.bpfreq   = [];      % to be filled by looping through foi_center
cfg_hfa.channel  = 'all';

pow = cell(size(fois));
maxs = zeros(size(fois));
for f = 1:numel(fois)
    cfg_hfa.bpfreq = bplim(f,:);%[fois(f_ix)-bws(f_ix)/2 fois(f_ix)+bws(f_ix)/2];
    fprintf('\n------> filtering: %.03f - %.03f\n', cfg_hfa.bpfreq(1), cfg_hfa.bpfreq(2));
    pow{f} = ft_preprocessing(cfg_hfa,roi);
    maxs(f) = max(pow{f}.trial{1});
end

%% Plot Raw Distributions
bins = linspace(0,max(maxs),500);
figure;
subplot(15,2,1);
histogram(roi.trial{1});
title(['orig: skew = ' num2str(skewness(roi.trial{1})) '; kurt = ' num2str(kurtosis(roi.trial{1}))]);

means = zeros(size(fois));
vars  = zeros(size(fois));
skews = zeros(size(fois));
kurts = zeros(size(fois));
for f = 1:numel(fois)
    means(f) = mean(pow{f}.trial{1});
    vars(f)  = std(pow{f}.trial{1});
    skews(f) = skewness(pow{f}.trial{1});
    kurts(f) = kurtosis(pow{f}.trial{1});
    if f>14; c=2; else; c=1;end
    sp = fn_rowcol2subplot_ix(15,2,mod(f,15)+1,c);
    subplot(15,2,sp);
    histogram(pow{f}.trial{1},bins);
    title([num2str(bplim(f,1)) '-' num2str(bplim(f,2)) ' Hz: skew = ' num2str(skews(f)) '; kurt = ' num2str(kurts(f))]);
end
% suptitle('Raw Power Distributions');

%% Plot Log Transforms
bins = linspace(-10,5,500);
figure;
subplot(15,2,1);
histogram(roi.trial{1});
title(['original: skew = ' num2str(skewness(roi.trial{1})) '; kurt = ' num2str(kurtosis(roi.trial{1}))]);

lmeans = zeros(size(fois));
lvars  = zeros(size(fois));
lskews = zeros(size(fois));
lkurts = zeros(size(fois));
for f = 1:numel(fois)
    lmeans(f) = mean(log(pow{f}.trial{1}));
    lvars(f)  = std(log(pow{f}.trial{1}));
    lskews(f) = skewness(log(pow{f}.trial{1}));
    lkurts(f) = kurtosis(log(pow{f}.trial{1}));
    if f>14; c=2; else; c=1;end
    sp = fn_rowcol2subplot_ix(15,2,mod(f,15)+1,c);
    subplot(15,2,sp);
    histogram(log(pow{f}.trial{1}),bins);
    title([num2str(bplim(f,1)) '-' num2str(bplim(f,2)) ' Hz: skew = ' num2str(skews(f)) '; kurt = ' num2str(kurts(f))]);
end

%%
fig_name = [SBJ '_log_power_' roi.label{1}];
figure('Name',fig_name);
subplot(2,2,1);hold on;
plot(fois,means,'r');
plot(fois,lmeans,'b');
xlabel('frequency');
ylabel('Mean');
legend('raw','log');
scatter(0,mean(roi.trial{1}),25,'k');
title('Mean');

subplot(2,2,2);hold on;
plot(fois,vars,'r');
plot(fois,lvars,'b');
xlabel('frequency');
ylabel('StD');
legend('raw','log');
%scatter(0,std(roi.trial{1}),25,'k');
title('Standard Deviation');

subplot(2,2,3);hold on;
plot(fois,skews,'r');
plot(fois,lskews,'b');
xlabel('frequency');
ylabel('Skew');
legend('raw','log');
scatter(0,skewness(roi.trial{1}),25,'k');
title('Skewness');

subplot(2,2,4);hold on;
plot(fois,kurts,'r');
plot(fois,lkurts,'b');
xlabel('frequency');
ylabel('Kurtosis');
legend('raw','log');
scatter(0,kurtosis(roi.trial{1}),25,'k');
title('Kurtosis');


end