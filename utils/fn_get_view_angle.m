function view_angle = fn_get_view_angle(hemi,roi_id)
%% Define the view angle given hemisphere and ROI
%   OFC = bottom

% If no ROI, assume lateral view
if isempty(roi_id)
    roi_id = 'LPFC';
end

if strcmp(roi_id,'OFC')
    view_angle = [0 -90];   % from the bottom
elseif strcmp(hemi,'l')
    if any(strcmp(roi_id,{'LPFC','INS','TMP','PAR','MTL','lat','deep','ROI'}))
        view_angle = [-90 0];    % from the left
    elseif strcmp(roi_id,'MPFC')
        view_angle = [90 0];    % from the right
    else
        error(['Bad combo of hemi (' hemi ') and roi_id (' roi_id ')']);
    end
elseif strcmp(hemi,'r')
    if any(strcmp(roi_id,{'LPFC','INS','TMP','PAR','MTL','lat','deep','ROI'}))
        view_angle = [90 0];    % from the right
    elseif strcmp(roi_id,'MPFC')
        view_angle = [-90 0];    % from the left
    else
        error(['Bad combo of hemi (' hemi ') and roi_id (' roi_id ')']);
    end
else
    error(['Bad combo of hemi (' hemi ') and roi_id (' roi_id ')']);
end

end