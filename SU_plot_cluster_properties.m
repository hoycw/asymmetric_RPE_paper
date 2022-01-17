if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(genpath([app_dir 'wave_clus/']));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Get filenames
raw_micro_dir = '/Volumes/hoycw_clust/PRJ_Error/data/IR75/00_raw/SU_2018-06-01_14-37-09/micro/';
cluster_dir = '/Volumes/hoycw_clust/PRJ_Error/data/IR75/02_preproc/micro_clusters/semi_auto/';
spike_dir = '/Volumes/hoycw_clust/PRJ_Error/data/IR75/02_preproc/micro_spikes/';
clust_fnames = dir([cluster_dir 'times*3.mat']);
spike_fnames = dir([spike_dir '*spikes.mat']);
!! rebuild this as a struct for each microwire!
 !! can use probe names from SBJ_vars!

%% Compute Stats
ch_name = cell(size(spike_fnames));
n_clust = zeros(size(spike_fnames));
n_clust_probe = zeros([6 1]);
probe_names = {'mlam','mlhh','mlth','mram','mrhh','mrth'};
probe_print = {'LAM','LHH','LTH','RAM','RHH','RTH'};
frate   = cell(size(clust_fnames));
n_spikes= cell(size(clust_fnames));
per_bad = cell(size(clust_fnames));
waves  = cell(size(clust_fnames));
for s = 1:numel(spike_fnames)
    uscores = strfind(spike_fnames(s).name,'_');
    ch_name{s} = spike_fnames(s).name(1:uscores(1)-1);
    if s==1
%         micro = ft_read_neuralynx_interp({[raw_micro_dir ch_name{s} '_0003.ncs']});
%         data_len = micro.time{1}(end);
    end
    
    % Check for clusters
    c = find(strcmp({clust_fnames.name},['times_' ch_name{s} '_0003.mat']));
    if ~isempty(c)
        load(clust_files(c).name);
        n_clust(s) = sum(unique(cluster_class(:,1))>0);
        n_clust_probe(strcmp(probe_names,ch_name{s}(1:end-1))) = n_clust_probe(strcmp(probe_names,ch_name{s}(1:end-1)))+n_clust(s);
        n_spikes{c} = zeros([n_clust(s) 1]);
        frate{c} = zeros([n_clust(s) 1]);
        waves{c} = zeros([n_clust(s) size(spikes,2)]);
        for c2 = 1:n_clust(s)
            n_spikes{c}(c2) = sum(cluster_class(:,1)==c2);
            frate{c}(c2) = n_spikes{c}(c2)/data_len;
            waves{c}(c2,:) = mean(spikes(cluster_class(:,1)==c2,:),1);
            spike_time_dif = diff(cluster_class(cluster_class(:,1)==c2,2));
            per_bad{c}(c2) = 100*sum(spike_time_dif<3)/sum(cluster_class(:,1)==c2);
        end
    end
end

%% Plot Stats
% Bar graph on yield per probe
figure;
bar(n_clust_probe);
set(gca,'XTickLabel',probe_print);
xlabel('Probe');
ylabel('# Units');
title(['Yield by Probe (total n = ' num2str(sum(n_clust_probe)) ')']);

% Plot histogram of contamination rates
figure;
all_per_bad = horzcat(per_bad{:});
bad_bins = linspace(0,4,20);
histogram(all_per_bad,bad_bins);
title('Contamination Rates across All Units');

% Plot histogram of firing rates
all_frates = vertcat(frate{:});
fr_bins = linspace(0,max(all_frates),20);
figure;
histogram(all_frates,fr_bins);
%violinplot(all_frates);
title('Firing Rate Distribution');

% Plot FR vs. Contamination
figure;
scatter(all_per_bad,all_frates);
xlabel('Contamination Rate');
ylabel('Firing Rate');

