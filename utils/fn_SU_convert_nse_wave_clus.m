function [] = fn_SU_convert_nse_wave_clus(SBJ)
error('did not check this runs since Stroop conversion, why use .nse anyways?');
%% Load NLX .nse spike file and prepare for wave_clus input
%   This script was developed based on the issue in the wave_clus github:
%           https://github.com/csn-le/wave_clus/issues/34

%% Paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%% Paths
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(genpath([app_dir 'wave_clus/']));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% SBJ vars
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);

%% Load data
% Find first micro to get time stamps
for b_ix = 1:numel(SBJ_vars.block_name)
    ncs_fnames = dir([SBJ_vars.dirs.nlx{b_ix} 'micro/*.ncs']);
    nse_fnames = dir([SBJ_vars.dirs.nlx{b_ix} 'spikes/*.nse']);
    hdr = ft_read_header({[SBJ_vars.dirs.nlx{b_ix} 'micro/' ncs_fnames(1).name]});
    
    for ch = 1:numel(nse_fnames)
        if nse_fnames(ch).bytes>17000   % skip if file is empty (16384 = no spikes)
            spikeImport = ft_read_spike([SBJ_vars.dirs.nlx 'spikes/' nse_fnames(ch).name]);
            spikes = squeeze(spikeImport.waveform{1}(1,:,:))';
            par = set_parameters;
            par.w_pre = spikeImport.hdr.AlignmentPt;      % This is the default for Cheetah and Pegasus, confirmed in IR67 WedPM_Pegasus2018-01-24_17-59-36.log
            par.w_post = size(spikes,2)-par.w_pre;
            par.sr = spikeImport.hdr.SamplingFrequency;
            par.channels = 1;
            
            Fs = spikeImport.hdr.SamplingFrequency;
            index = 1000*double((spikeImport.timestamp{1}-spikeImport.hdr.FirstTimeStamp))/hdr.TimeStampPerSample/hdr.Fs;
            
            out_dir = [SBJ_vars.dirs.import 'mat_spikes/'];
            if ~exist(out_dir,'dir'); mkdir(out_dir); end
            save([out_dir nse_fnames(ch).name(1:end-4) '.mat'],'-v7.3','spikes','par','index');
            if ~exist([out_dir 'nse_wave_clus_param.mat'])
                save([out_dir 'nse_wave_clus_param.mat'],'-v7.3','par');
            end
        end
    end
end

end