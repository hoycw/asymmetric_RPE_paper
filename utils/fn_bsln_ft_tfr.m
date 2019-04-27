function bslnd_tfr = fn_bsln_ft_tfr(tfr, bsln_lim, bsln_type, n_boots)
%% Baseline correct one TFR based on bsln_lim epoch, both from ft_freqanalysis
% INPUTS:
%   tfr [ft dataset] - full output of ft_freqanalysis
%   bsln_lim [int, int]- 2 int array of TIME indices for [start, end] of baseline period
%   bsln_type [str]   - type of baseline to implement
%       'zboot'  = pool all baselines, bootstrap means+SDs, z-score all
%           trials to the mean+SD of the bootstrapped distribution
%       'zscore' = subtract mean and divide by SD
%       'demean' = subtract mean
%       'my_relchange' = subtract mean, divide by mean (results in % change)
% OUTPUTS:
%   bslnd_tfr [ft dataset] - same tfr but baseline corrected
[~, app_dir] = fn_get_root_dir();
addpath([app_dir 'fieldtrip/']);
ft_defaults
rng('shuffle'); % seed randi with time

if ~strcmp(tfr.dimord,'rpt_chan_freq_time')
    error('Check dimord to be sure trial dimension is first!')
end

% Select baseline data
cfgs = [];
cfgs.latency = bsln_lim;
bsln_tfr = ft_selectdata(cfgs,tfr);

bslnd_tfr = tfr;
for ch = 1:size(tfr.powspctrm,2)
    for f = 1:size(tfr.powspctrm,3)
        % Create bootstrap distribution if necessary
        if strcmp(bsln_type,'zboot')
            sample_means = NaN([1 n_boots]);
            sample_stds  = NaN([1 n_boots]);
            for boot_ix = 1:n_boots
                % Select a random set of trials (sampling WITH REPLACEMENT!)
                shuffle_ix = randi(size(tfr.powspctrm,1),[1 size(tfr.powspctrm,1)]);
                % Select baseline data and compute stats
                bsln_data = bsln_tfr.powspctrm(shuffle_ix,ch,f,:);
                sample_means(boot_ix) = nanmean(bsln_data(:));
                sample_stds(boot_ix)  = nanstd(bsln_data(:));
            end
        end
        
        % Perform Baseline Correction
        for t = 1:size(tfr.powspctrm,1)
            trials  = tfr.powspctrm(t,ch,f,:);
            trl_bsln    = bsln_tfr.powspctrm(t,ch,f,:);
            switch bsln_type
                case 'zboot'
                    bslnd_tfr.powspctrm(t,ch,f,:) = (trials-mean(sample_means))/mean(sample_stds);                    
                case 'zscore'
                    bslnd_tfr.powspctrm(t,ch,f,:) = (trials-nanmean(trl_bsln))/nanstd(trl_bsln);
                case 'demean'
                    bslnd_tfr.powspctrm(t,ch,f,:) = trials-nanmean(trl_bsln);
                case 'my_relchange'
                    bslnd_tfr.powspctrm(t,ch,f,:) = (trials-nanmean(trl_bsln))/nanmean(trl_bsln);
                otherwise
                    error(['Unknown bsln_type: ' bsln_type])
            end
        end
    end
end

end
