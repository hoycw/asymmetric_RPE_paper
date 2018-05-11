function [labels, colors, einfo_col] = fn_roi_label_styles(roi_id)
%% Converts the name of a set of ROIs into labels, plotting colors
% condition_name: [str] 'ROI', 'gROI', 'INS', 'LPFC', 'MPFC', 'OFC'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%
% einfo_col = 2 for specific ROIs, 3 for general ROIs

% if length(cond_lab) == 1
switch roi_id
    case 'ROI'
        load('~/PRJ_Error/data/full_roi_lists.mat');
        labels = all_rois;
        % Exclude FWM, '', OUT
        labels(strmatch('FWM',labels,'exact')) = [];
        labels(strmatch('OUT',labels,'exact')) = [];
        labels(strmatch('',labels,'exact')) = [];
        einfo_col = 2;
    case 'gROI'
        labels = {'LPFC','MPFC','INS','OFC'};
        einfo_col = 3;
    case 'thryROI'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA','daINS','vaINS','FO'};
        einfo_col = 2;
    case 'LPFC'
        labels = {'FPC','DLPFC','VLPFC','PM','M1','S1'};
        einfo_col = 2;
    case 'MPFC'
        labels = {'ACC','preSMA','aMCC','SMA','pMCC','Precuneus'};
        einfo_col = 2;
    case 'INS'
        labels = {'vaINS','daINS','FO','mINS','pINS'};
        einfo_col = 2;
    case 'OFC'
        labels = {'mOFC','lOFC'};
        einfo_col = 2;
    case 'ALL'
        load('~/PRJ_Error/data/full_roi_lists.mat');
        labels = all_rois;
        einfo_col = 2;
    otherwise
        error(strcat('Unknown roi_id: ',roi_id));
end

% Get colors
colors = cell(size(labels));
for roi_ix = 1:numel(labels)
    colors{roi_ix} = fn_roi2color(labels{roi_ix});
end

end
