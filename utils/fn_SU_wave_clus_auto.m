function fn_SU_wave_clus_auto(SBJ,block_ix,varargin)
%% Run wave_clus automated spike sorting
% Paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end
addpath(genpath([app_dir 'wave_clus/']));

%% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'chan_lab')
            chan_lab = varargin{v+1};
        elseif strcmp(varargin{v},'spike_mat_dir')
            spike_mat_dir = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

if ~exist('chan_lab','var')
    chan_lab = 'all';
end

%% Do clustering
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
if ~exist('spike_mat_dir','var') % .ncs files (default)
    cd([SBJ_vars.dirs.nlx{block_ix} 'micro/']);
    Get_spikes(chan_lab);
    Do_clustering(chan_lab);
else            % .mat files
    cd(spike_mat_dir);
    load('nse_wave_clus_param.mat');
    % Get list of spike files
    mat_fnames = dir([spike_mat_dir '/*mat']);
    mat_fnames(strcmp({mat_fnames.name},'nse_wave_clus_param.mat')) = [];
    Get_spikes({mat_fnames.name},'par',par);
    Do_clustering(chan_lab);
end

end