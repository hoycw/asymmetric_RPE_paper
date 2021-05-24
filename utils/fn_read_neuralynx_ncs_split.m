function [ncs_out] = fn_read_neuralynx_ncs_split(fname)
%% Read neuralynx .ncs file using Fieldtrip read_neuralynx_ncs
% Based on ft_read_neuralynx_interp from Martink Vinck and Fieldtrip
error('dont use this, use fn_read_neuralynx_interp');
if exist('/home/knight/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Check header info
% first check if these are indeed ncs files
ftype = zeros(length(fname), 1);
for i=1:length(fname)
  if     ft_filetype(fname{i}, 'neuralynx_ncs')
    ftype(i) = 1;
  end
end
if ~all(ftype==1), error('some files do not correspond to ncs files'); end

% number of channels
nchans = length(fname);

% get the original headers
for i=1:nchans
  orig(i) = ft_read_header(fname{i});
end

% check if they all have the same sampling frequency, otherwise return error
for i=1:length(orig), SamplingFrequency(i) = orig(i).Fs; end
if any(SamplingFrequency~=SamplingFrequency(1))
  error('not all channels have the same sampling rate');
end

%% check for invalid samples
bad_rec_idx = cell([nchans 1]);
for i=1:nchans
    ncs = fn_read_neuralynx_ncs(fname{i});
    % count invalid samples per channel
    bad_rec_idx{i} = find(ncs.NumValidSamp(1:end-1)<512);
    % check for alignment across channels
    if i>1 && ~all(bad_rec_idx{i}==bad_rec_idx{1})
        error('not all channels have same records with invalid samples');
    end
end

% Find which are true gaps vs. repeats
gap_size = diff(ncs.TimeStamp);
gap_sizes = unique(gap_size);
good_gap_len  = [];
for g_ix = 1:numel(gap_sizes)
    if sum(gap_size==gap_sizes(g_ix))>1
        good_gap_len = [good_gap_len gap_sizes(g_ix)];
    end
end
% find bad gaps less then godo gap size
repeat_gap_idx  = [];
long_gap_idx = [];
for g_ix = 1:numel(gap_sizes)
    if ~any(gap_sizes(g_ix)==good_gap_len)
        if all(gap_sizes(g_ix)<good_gap_len)
            repeat_gap_idx = [repeat_gap_idx find(gap_size==gap_sizes(g_ix))];
        else
            long_gap_idx = [long_gap_idx find(gap_size==gap_sizes(g_ix))];
        end
    end
end
if ~isempty(setdiff(bad_rec_idx{1},[long_gap_idx repeat_gap_idx]))
    error('check invalid samples!');
end
record_ts_gap  = (ncs.TimeStamp(2)-ncs.TimeStamp(1))/512;
repeat_timestamps = ncs.TimeStamp(repeat_gap_idx);
long_gap_timestamps = ncs.TimeStamp(long_gap_idx);

%% Load data
% ncs_out = read_neuralynx_ncs(filename);
for i = 1:nchans
    cfg         = [];
    cfg.dataset = fname{i};
    data        = ft_preprocessing(cfg);
    ts          = ft_read_data(cfg.dataset, 'timestamp', true);
    ts_orig     = ts;
    
    long_gap_ts_ix = zeros(size(long_gap_timestamps));
    for g_ix = 1:numel(long_gap_timestamps)
        long_gap_ts_ix(g_ix) = find(ts_orig==long_gap_timestamps(g_ix));
    end
    repeat_gap_ts_ix = cell(size(repeat_timestamps));
    for g_ix = 1:numel(repeat_timestamps)
        repeat_gap_ts_ix{g_ix} = find(ts_orig==repeat_timestamps(g_ix));
    end
    
    % original timestamaps in doubles, with the minimum ts subtracted
    ts = double(ts_orig-ts_orig(1));
    ts_step = ts(2)-ts(1);
    % convert from us to s
    if ts_step*10^-6~=1/SamplingFrequency; error('sampling rate mismatch!'); end
    ts = ts*10^-6;
    
    % replace time with better timestamps
    data.time{1} = ts;
    long_gap_times = ts(long_gap_ts_ix);
    
    % Split data
    chunks = cell([numel(long_gap_times)+1 1]);
    cfgs   = [];
    cfgs.latency = [ts(1) long_gap_times(1)];
    chunks{1}    = ft_selectdata(cfgs,data);
    for chunk_ix = 2:numel(long_gap_times)
        cfgs.latency = [long_gap_times(chunk_ix-1)+ts_step long_gap_times(chunk_ix)];
        chunks{chunk_ix} = ft_selectdata(cfgs,data);
        fprintf('in loop ix = %d\n',chunk_ix);
    end
    cfgs.latency = [long_gap_times(end)+ts_step ts(end)];
    chunks{numel(chunks)} = ft_selectdata(cfgs,data);
    
end

end
