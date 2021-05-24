function good_lab = fn_select_elec_lab_match(elec, hemi, atlas_id, roi_id)
%% Apply selection criteria and return list of matching elec labels
% INPUTS:
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
% OUTPUTS:
%   good_lab [cell array strs] - list of labels matching all criteria

% Start with all labels
lab_match = true(size(elec.label));

% Match hemisphere
if ~strcmp(hemi,'b')
    lab_match = all([lab_match strcmp(elec.hemi,hemi)],2);
end

% Match Atlas
if ~isempty(atlas_id)
    if ~isempty(roi_id)
        [roi_list, ~, roi_field] = fn_roi_label_styles(roi_id);
        roi_match = zeros([numel(elec.label) numel(roi_list)]);
        for roi_ix = 1:numel(roi_list)
            roi_match(:,roi_ix) = strcmp(elec.(roi_field),roi_list(roi_ix));
        end
        lab_match = all([lab_match any(roi_match,2)],2);
    else
        lab_match = all([lab_match ~strcmp(elec.atlas_lab,'no_label_found')],2);
    end
end

if ~all(lab_match)
    fprintf('\themi (%s) and roi_id (%s): tossed %i/%i elecs\n',...
               hemi,roi_id,sum(~lab_match),numel(elec.label));
else
    fprintf('hemi (%s) and roi_id (%s): keeping all elecs\n',hemi,roi_id);
end

good_lab = elec.label(lab_match);

end
