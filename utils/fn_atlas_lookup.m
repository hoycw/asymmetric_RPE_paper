function [elec_lab] = fn_atlas_lookup(elec, atlas, varargin)
% generate a list of roi_labels for a given elec struct and atlas
% INPUTS:
%   elec [struct] - elec struct
%   atlas [struct] - atlas loaded from ft_read_atlas
%   varargin [cell array] - used to pass in {'max_qry_rng','min_qry_range'}
%       either should be pased in form 'm**_qry_rng',[1 3 5]


% edit generate_electable.m % more extensive version

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'max_qry_rng') && any(ismember(varargin{v+1},[1 3 5]))
            max_qry_rng = varargin{v+1};
        elseif strcmp(varargin{v},'min_qry_rng') && any(ismember(varargin{v+1},[1 3 5]))
            min_qry_rng = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

if ~exist('max_qry_rng','var')
    max_qry_rng = 5;
end
if ~exist('min_qry_rng','var')
    min_qry_rng = 1;
end

%% Find ROI Labels of electrodes (Create elec table)
% Initialize the new fields
elec_lab = elec;
elec_lab.atlas_id     = atlas.name;
elec_lab.atlas_lab    = cell(size(elec.label));
elec_lab.atlas_prob   = zeros(size(elec.label));
elec_lab.atlas_qryrng = zeros(size(elec.label));
elec_lab.atlas_qrysz  = zeros(size(elec.label));
elec_lab.atlas_lab2   = cell(size(elec.label));
elec_lab.atlas_prob2  = cell(size(elec.label));

cfg = [];
cfg.atlas              = atlas;
cfg.inputcoord         = elec.coordsys;
cfg.output             = 'multiple';
cfg.minqueryrange      = min_qry_rng;
cfg.maxqueryrange      = max_qry_rng;
cfg.querymethod        = 'sphere';
cfg.round2nearestvoxel = 'yes';
for e = 1:numel(elec.label)
    % Search for this elec
    cfg.roi = elec.chanpos(e,:);
    report = ft_volumelookup(cfg,atlas);
    
    % Report label for this elec
    [sorted_cnt, cnt_idx] = sort(report.count,1,'descend');
    match_cnt = find(sorted_cnt);
    
    % Assign highest match as label
    if numel(match_cnt)>=1
        elec_lab.atlas_lab{e}  = report.name{cnt_idx(1)};
        elec_lab.atlas_prob(e) = report.count(cnt_idx(1))/sum(report.count);
        % report.usedqueryrange bug work around:
        %   ft_volumelookup only assigns usedQR to last ROI in list
        %   if best match isn't last, usedQR is empty --> error assigning to elec_lab
        %   instead: grab whatever the usedQR for that iteration was (should be the same)
        if numel(report.usedqueryrange{cnt_idx(1)})~=1
            if numel([report.usedqueryrange{:}])>1      % check that there was only one usedQR
                error('work around for missing usedQR didnt work!');
            end
            elec_lab.atlas_qryrng(e) = [report.usedqueryrange{:}];
        else
            elec_lab.atlas_qryrng(e) = report.usedqueryrange{cnt_idx(1)};
        end
        elec_lab.atlas_qrysz(e) = sum(report.count);
        % add additional labels if needed
        if numel(match_cnt)>1
            elec_lab.atlas_lab2{e}  = report.name(cnt_idx(2:numel(match_cnt)));
            elec_lab.atlas_prob2{e} = report.count(cnt_idx(2:numel(match_cnt)))/sum(report.count);
        else
            elec_lab.atlas_lab2{e}  = '';
            elec_lab.atlas_prob2{e} = NaN;
        end
    else   % No matches found
        error(['No matches found for elec: ' elec.label{e}]);
    end
%     for j=1:length(match_cnt);
%         found{j} = report.name{cnt_idx(j)};%[num2str(report.count(ind(j))) ': ' report.name{ind(j)}];
%     end
%     roi_labels{e} = strjoin(found,',');
end

% nf_cnt = 0;
% for e = 1:numel(elec.label)
%     disp([num2str(e) ': ' elec.label{e} ' = ' roi_labels{e}]);
%     if ~isempty(strfind(roi_labels{e},'no_label_found'))
%         nf_cnt= nf_cnt+1;
%     end
% end
% disp(['nf_cnt = ' num2str(nf_cnt)]);
% % Import atlas of interest
% atlas = ft_read_atlas([ft_dir 'template/atlas/aal/ROI_MNI_V4.nii']);
% [~, indx] = max(labels.count);
% labels.name(indx)

%% Check that all atlas_prob add to 1
for e = 1:numel(elec_lab.label)
    if elec_lab.atlas_prob(e)+sum(elec_lab.atlas_prob2{e})<0.99999  % sometimes it's 0.99999999999999988898
        error(['Electrode ' elec_lab.label{e} ' has atlas_prob = '...
            num2str(elec_lab.atlas_prob(e)+sum(elec_lab.atlas_prob2{e}))]);
    end
end

end


