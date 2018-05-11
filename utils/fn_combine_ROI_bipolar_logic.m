function [new_roi, new_groi, new_roi2, new_tissue, gm_weight, out] = fn_combine_ROI_bipolar_logic(...
                            roi1,groi1,tiss_p1,tiss_s1,roi2,groi2,tiss_p2,tiss_s2)
% Logic to determine the properties assigned to a bipolar referenced pair of electrodes
% INPUTS: two versions of each for the two electrodes in the pair
%   roi - detailed roi
%   groi - gross level roi (LPFC, MPFC, OFC, INS, FWM=frontal white matter)
%   tiss_p - primary tissue type (GM, WM, OUT)
%   tiss_s - secondary tissue type
% OUTPUTS:
%   new_roi - specific ROI for this pair
%   new_groi - general ROI for this pair
%   new_roi2 - secondary ROI contributing to this pair
%   new_tissue - primary tissue type contributing to this pair
%   gm_weight - portion of this pair with GM signal
%   out - 0/1 binary flag if there's any out of brain signal in this pair

%% Tissue Logic
gm_weight1 = fn_gm_weight(tiss_p1,tiss_s1);
gm_weight2 = fn_gm_weight(tiss_p2,tiss_s2);
gm_weight  = mean([gm_weight1 gm_weight2]);
if gm_weight >= 0.5
    new_tissue = 'GM';
else
    new_tissue = 'WM';
end

%% ROI Logic
if strcmp(roi1,roi2)            % Both are the same ROI
    roi_ix = 1;
else
    if gm_weight1 > gm_weight2  % One is more in GM than the other
        roi_ix = 1;
    elseif gm_weight2 > gm_weight1
        roi_ix = 2;
    else                        % Both equally in GM
        if strcmp(roi1,'FWM')   % One is in deep white matter
            roi_ix = 2;
        else                    % Either roi2 is FWM or it's even, so just choose deepest contact
            roi_ix = 1;
        end
    end
end

%% Select Info from Correct ROI
if roi_ix == 1
    new_roi  = roi1;
    new_groi = groi1;
    new_roi2 = roi2;
else
    new_roi  = roi2;
    new_groi = groi2;
    new_roi2 = roi1;
end

%% Flag if any potential out of brain signal
if any(strcmp({tiss_p1,tiss_s1,tiss_p2,tiss_s2},'OUT'))
    out = 1;
else
    out = 0;
end

end


function [gm_weight] = fn_gm_weight(primary_tissue,secondary_tissue)
%% Computes the portion of a bipolar electrode pair with GM signal
if strcmp(primary_tissue,'GM')
    if strcmp(secondary_tissue,'')  % GM, GM
        gm_weight = 1;
    else                            % GM, WM; GM, OUT
        gm_weight = 0.6;
    end
else
    if strcmp(secondary_tissue,'')  % WM, WM
        gm_weight = 0;
    else                            % WM, GM; OUT, GM
        gm_weight = 0.4;
    end
end
end