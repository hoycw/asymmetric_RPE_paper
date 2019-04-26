function db_out = SBJ00b_view_preclean(SBJ,block_ix,keep_db_out,varargin)
%% View SBJ00_preclean output with channels colored by type
% INPUTS:
%   SBJ [str] - name of patient
%   block_ix [int] - which block to plot
%   keep_db_out [0/1] - return the output of ft_databrowser to save epochs
%   varargin:
%       bad_epochs [str or Nx2 int] - 'load'; or [start stop] indices
%       reorder [cell array] - either list of labels or {} for alphabetical
%       ylim [1x2 int] - [min max] for ft_databrowser ylims
% OUTPUTS:
%   db_out [struct] - output of ft_databrowser if keep_db_out==1

%% Assign colors
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain
%   all the rest = 4
bad_color = [1 0 0];
sus_color = [0.7 0 0.7];
out_color = [0 0.1 1];
all_color = [0.3 0.3 0.4];
colormap  = [bad_color; sus_color; out_color; all_color];

%% Process varargin
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'bad_epochs')
            if ischar(varargin{v+1}) || size(varargin{v+1},2)==2
                bad_epochs = varargin{v+1};
            else
                error('unexpected bad_epochs input');
            end
        elseif strcmp(varargin{v},'reorder') && iscell(varargin{v+1})
            new_lab_order = varargin{v+1};
        elseif strcmp(varargin{v},'ylim') && numel(varargin{v+1})==2
            y_lim = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'hoycw/Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% ========================================================================
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% %% Keep any real but "bad" channels
% toss = ones(size(SBJ_vars.ch_lab.bad));
% for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
%     ch_match = zeros(size(SBJ_vars.ch_lab.probes));
%     for probe_ix = 1:numel(SBJ_vars.ch_lab.probes)
%         if ~isempty(strfind(SBJ_vars.ch_lab.bad{bad_ix},SBJ_vars.ch_lab.probes{probe_ix}))
%             ch_match(probe_ix) = 1;
%         end
%     end
%     if any(ch_match)
%         toss(bad_ix) = 0;
%     end
% end
% SBJ_vars.ch_lab.bad = SBJ_vars.ch_lab.bad(logical(toss));

%% Process channel labels
% Handle prefix and suffix
if isfield(SBJ_vars.ch_lab,'prefix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)        
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix}];
    end
    for ref_ix = 1:numel(SBJ_vars.ch_lab.ref_exclude)
        SBJ_vars.ch_lab.ref_exclude{ref_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ref_exclude{ref_ix}];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eeg{eeg_ix}];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eog{eog_ix}];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.photod{1}]};
end
if isfield(SBJ_vars.ch_lab,'suffix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix];
    end
    for ref_ix = 1:numel(SBJ_vars.ch_lab.ref_exclude)
        SBJ_vars.ch_lab.ref_exclude{ref_ix} = [SBJ_vars.ch_lab.ref_exclude{ref_ix} SBJ_vars.ch_lab.suffix];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.eeg{eeg_ix} SBJ_vars.ch_lab.suffix];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.eog{eog_ix} SBJ_vars.ch_lab.suffix];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.suffix]};
end
% bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);
% eeg_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eeg);
% eog_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eog);
photod_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.photod);

%% Load data
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{block_ix});
else
    block_suffix = '';
end

load([SBJ_vars.dirs.preproc SBJ '_preclean' block_suffix '.mat']);

%% Select Neural Data
junk_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad(SBJ_vars.ch_lab.bad_code==0));

cfg = [];
cfg.channel = [{'all','-EDF Annotations'},junk_ch_neg,photod_ch_neg];
data = ft_selectdata(cfg,data);

% Name EEG/EOG to stick together, move to bottom of sort, add # if necessary for sorting
for l = 1:numel(data.label)
    if any(strcmp(data.label{l},SBJ_vars.ch_lab.eeg))
        data.label{l} = ['zEEG0-' data.label{l}];
    elseif any(strcmp(data.label{l},SBJ_vars.ch_lab.eog))
        data.label{l} = ['zEOG0-' data.label{l}];
    end
end
for l = 1:numel(SBJ_vars.ch_lab.eeg)
    SBJ_vars.ch_lab.eeg{l} = ['zEEG0-' SBJ_vars.ch_lab.eeg{l}];
end
for l = 1:numel(SBJ_vars.ch_lab.eog)
    SBJ_vars.ch_lab.eog{l} = ['zEOG0-' SBJ_vars.ch_lab.eog{l}];
end

%% Reorder
if exist('new_lab_order','var')
    data = fn_reorder_data(data, new_lab_order);
end

%% Plotting Set up
load([root_dir 'PRJ_Error/scripts/utils/cfg_plot.mat']);
if exist('bad_epochs','var')
    if ischar(bad_epochs) && strcmp(bad_epochs,'load')
        load([SBJ_vars.dirs.events SBJ '_bad_epochs_preclean' block_suffix '.mat']);
    end
    cfg_plot.artfctdef.visual.artifact = bad_epochs;
end
if exist('y_lim','var')
    cfg_plot.ylim = y_lim;
end

% Set the colors
cfg_plot.channelcolormap = colormap;

% Assign channels to colorgroups
cgroup = zeros([numel(data.label) 1]);
for l = 1:numel(data.label)
    if any(strcmp(data.label{l},SBJ_vars.ch_lab.bad))
        cgroup(l)  = SBJ_vars.ch_lab.bad_code(strcmp(data.label{l},SBJ_vars.ch_lab.bad));
    elseif any(strcmp(data.label{l},SBJ_vars.ch_lab.ref_exclude))
        cgroup(l)  = find(strcmp(SBJ_vars.ch_lab.bad_type,'sus'));
    elseif any(strcmp(data.label{l},SBJ_vars.ch_lab.eeg)) || any(strcmp(data.label{l},SBJ_vars.ch_lab.eog))
        cgroup(l)  = find(strcmp(SBJ_vars.ch_lab.bad_type,'out'));
    else
        cgroup(l)  = 4;
    end
    if cgroup(l)==4
        fprintf('%s \t all\n',data.label{l});
    else
        fprintf(2, '%s \t %s\n',data.label{l},SBJ_vars.ch_lab.bad_type{cgroup(l)});
    end
end
cfg_plot.colorgroups = cgroup;

%% Plot the data
if keep_db_out
    db_out = ft_databrowser(cfg_plot,data);
else
    ft_databrowser(cfg_plot,data);
end

%% Save
% save([SBJ_vars.dirs.events SBJ '_bad_epochs_preclean.mat'],'-v7.3','bad_epochs');

end
