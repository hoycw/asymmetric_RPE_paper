function [RGB] = fn_roi2color(roi)
%% Returns the RGB color code to plot a given ROI
% INPUTS:
%   roi [str or cell] - name (or list of names) of the ROI
%       can be specific or general (e.g., DLPFC or LPFC)
% OUTPUTS:
%   RGB [3 floats] - RGB for this ROI
%
% Values taken from: https://www.w3schools.com/colors/colors_groups.asp

if ischar(roi)
    n_rois = 1;
    roi = {roi};
elseif iscell(roi)
    n_rois = numel(roi);
end

RGB = zeros([n_rois 3]);
for r = 1:n_rois
    switch roi{r}
        % General ROIs
        case 'LPFC'         % Orange
            RGB(r,:) = [253 192 134]./256;
        case 'MPFC'         % Purple
            RGB(r,:) = [190 174 212]./256;
        case 'INS'          % Green
            RGB(r,:) = [127 201 127]./256;
        case 'OFC'          % Dark Blue
            RGB(r,:) = [56 108 176]./256;
        case 'PAR'          % Maroon
            RGB(r,:) = [128 0 0]./256;
        case 'TMP'          % Yellow
            RGB(r,:) = [1 1 0];%[0.2 0.2 0.2];
        case 'MTL'          % Cyan
            RGB(r,:) = [0 1 1];
        case 'OCC'          % Dark Gray
            RGB(r,:) = [0.2 0.2 0.2];
        case 'FWM'          % White
            RGB(r,:) = [1 1 1];
        case 'TWM'          % White
            RGB(r,:) = [1 1 1];
            
            % LPFC Subregions - reds
        case 'FPC'
            RGB(r,:) = [179 0 0]./256;
        case 'DLPFC'
            RGB(r,:) = [227 74 51]./256;
        case 'VLPFC'
            RGB(r,:) = [252 141 89]./256;
        case 'PM'
            RGB(r,:) = [253 187 132]./256;
        case 'M1'
            RGB(r,:) = [253 212 158]./256;
            %    case 'S1'
            %        RGB(r,:) = [254 240 217]./256;
            
            % MPFC Subregions - blues
        case 'ACC'
            RGB(r,:) = [8 81 156]./256;
            RGB(r,:) = [190 174 212]./256;
        case 'preSMA'
            RGB(r,:) = [49 130 189]./256;
            RGB(r,:) = [190 174 212]./256;
        case 'aMCC'
            RGB(r,:) = [107 174 214]./256;
            RGB(r,:) = [190 174 212]./256;
        case 'SMA'
            RGB(r,:) = [158 202 225]./256;
            RGB(r,:) = [190 174 212]./256;
        case 'pMCC'
            RGB(r,:) = [198 219 239]./256;
            RGB(r,:) = [190 174 212]./256;
        case 'PCC'
            RGB(r,:) = [239 243 255]./256;
            RGB(r,:) = [190 174 212]./256;
            
            % Insula Subregions - greens
        case 'vaINS'
            RGB(r,:) = [0 109 44]./256;
        case 'daINS'
            RGB(r,:) = [49 163 84]./256;
        case 'FO'
            RGB(r,:) = [116 196 118]./256;
        case 'mINS'
            RGB(r,:) = [161 217 155]./256;
        case 'pINS'
            RGB(r,:) = [199 233 192]./256;
            
            % OFC Subregions - yellows
        case 'mOFC'
            RGB(r,:) = [56 108 176]./256;
        case 'lOFC'
            RGB(r,:) = [56 108 176]./256;
            
            % Parietal - purples
        case 'IPL'
            RGB(r,:) = [129 15 124]./256;
        case 'SPL'
            RGB(r,:) = [136 86 167]./256;
        case 'S1'
            RGB(r,:) = [140 150 198]./256;
        case 'Precuneus'
            RGB(r,:) = [179 205 227]./256;
            
            % MTL (Medial Temporal Lobe) - bright colors (for now?)
        case 'HPC'
            RGB(r,:) = [0 1 1];
        case 'AMG'
            RGB(r,:) = [1 0 1];
        case 'STS'
            RGB(r,:) = [1 1 0];%[0.2 0.2 0.2];

            % Weird Cases
        case ''
            RGB(r,:) = [0 0 0];
        case 'OUT'          % Black
            warning('WARNING: Why are you trying to plot data mainly out of the brain???');
            RGB(r,:) = [0 0 0];
            
%             % Yeo Atlases
        case 'Vis'
            RGB(r,:) = [120  18 134]./256;
        case 'SM'
            RGB(r,:) = [70 130 180]./256;
        case 'DAttn'
            RGB(r,:) = [0 118  14]./256;
        case 'VAttn'
            RGB(r,:) = [196  58 250]./256;
        case 'Limb'
            RGB(r,:) = [220 248 164]./256;
        case 'FP'
            RGB(r,:) = [230 148  34]./256;
        case 'Def'
            RGB(r,:) = [205  62  78]./256;

            % Tissue Types
        case 'GM'
            RGB(r,:) = [0.4 0.4 0.4];
        case 'WM'
            RGB(r,:) = [1 1 1];
        case 'CSF'
            RGB(r,:) = [0 0 0];
        % case 'OUT' is covered above
    end
end

end
