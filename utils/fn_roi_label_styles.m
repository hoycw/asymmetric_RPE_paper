function [labels, colors, field] = fn_roi_label_styles(roi_id)
%% Converts the name of a set of ROIs into labels, plotting colors
% condition_name: [str] 'ROI', 'gROI', 'INS', 'LPFC', 'MPFC', 'OFC'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%
% einfo_col = 2 for specific ROIs, 3 for general ROIs

switch roi_id
    case 'ROI'
        load('~/PRJ_Error/data/full_roi_lists.mat');
        labels = all_rois;
        % Exclude FWM, '', OUT
        labels(strmatch('FWM',labels,'exact')) = [];
        labels(strmatch('TWM',labels,'exact')) = [];
        labels(strmatch('OUT',labels,'exact')) = [];
        labels(strmatch('',labels,'exact')) = [];
        field = 'ROI';
    case 'Yeo7'
        labels = {'Vis','SM','DAttn','VAttn','Limb','FP','Def'};
    case 'MPFCINS'
        labels = {'MPFC','INS'};
        field = 'gROI';
    case 'main3'
        labels = {'LPFC','MPFC','INS'};
        field = 'gROI';
    case 'mgROI'
        labels = {'LPFC','MPFC','INS','OFC','SM'};
        field = 'gROI';
    case 'gROI'
        labels = {'LPFC','MPFC','INS','OFC','SM','PAR','TMP','AMG','HPC'};%,'OCC'};
        field = 'gROI';
    case 'lat'
        labels = {'LPFC','PAR','TMP','OCC'};
        field = 'gROI';
    case 'deep'
        labels = {'INS','HPC','AMG'};
        field = 'gROI';
    case 'mnLPFC'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA'};
        field = 'ROI';
    case 'thryROI'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA','daINS','vaINS','FO'};
        field = 'ROI';
    case 'PAR'
        labels = {'S1','SPL','IPL','Precuneus'};
        field = 'ROI';
    case 'TMP'
        labels = {'STS'};
        field = 'ROI';
    case 'LPFC'
        labels = {'FPC','DLPFC','VLPFC','PM','M1'};
        field = 'ROI';
    case 'MPFC'
        labels = {'ACC','preSMA','aMCC','SMA','pMCC'};
        field = 'ROI';
    case 'INS'
        labels = {'vaINS','daINS','FO','mINS','pINS'};
        field = 'ROI';
    case 'OFC'
        labels = {'mOFC','lOFC'};
        field = 'ROI';
    case 'MTL'
        labels = {'AMG','HPC'};
        field = 'ROI';
    case {'tissue', 'tissueC'}
        labels = {'GM','WM','CSF','OUT'};
        field = 'tissue';
    case 'all'
        load('~/PRJ_Error/data/full_roi_lists.mat');
        labels = all_rois;
        field = 'ROI';
    otherwise
        error(strcat('Unknown roi_id: ',roi_id));
end

% Get colors
colors = cell(size(labels));
for roi_ix = 1:numel(labels)
    colors{roi_ix} = fn_roi2color(labels{roi_ix});
end

end
