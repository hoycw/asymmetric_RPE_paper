function [labels, colors] = fn_roi_label_styles(roi_id)
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
        labels(strmatch('TWM',labels,'exact')) = [];
        labels(strmatch('OUT',labels,'exact')) = [];
        labels(strmatch('',labels,'exact')) = [];
    case 'Yeo7'
        labels = {'Vis','SM','DAttn','VAttn','Limb','FP','Def'};
    case 'main3'
        labels = {'LPFC','MPFC','INS'};
    case 'mgROI'
        labels = {'LPFC','MPFC','INS','OFC'};
    case 'gROI'
        labels = {'LPFC','MPFC','INS','OFC','PAR','TMP','AMG','HPC'};%,'OCC'};
    case 'lat'
        labels = {'LPFC','PAR','TMP','OCC'};
    case 'deep'
        labels = {'INS','HPC','AMG'};
    case 'mnLPFC'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA'};
    case 'thryROI'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA','daINS','vaINS','FO'};
    case 'PAR'
        labels = {'S1','SPL','IPL','Precuneus'};
    case 'TMP'
        labels = {'STS'};
    case 'LPFC'
        labels = {'FPC','DLPFC','VLPFC','PM','M1'};
    case 'MPFC'
        labels = {'ACC','preSMA','aMCC','SMA','pMCC'};
    case 'INS'
        labels = {'vaINS','daINS','FO','mINS','pINS'};
    case 'OFC'
        labels = {'mOFC','lOFC'};
    case 'MTL'
        labels = {'AMG','HPC'};
    case {'tissue', 'tissueC'}
        labels = {'GM','WM','CSF','OUT'};
    case 'all'
        load('~/PRJ_Error/data/full_roi_lists.mat');
        labels = all_rois;
    otherwise
        error(strcat('Unknown roi_id: ',roi_id));
end

% Get colors
colors = cell(size(labels));
for roi_ix = 1:numel(labels)
    colors{roi_ix} = fn_roi2color(labels{roi_ix});
end

end
