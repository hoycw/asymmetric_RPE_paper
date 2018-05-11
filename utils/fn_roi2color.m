function [RGB] = fn_roi2color(roi)
%% Returns the RGB color code to plot a given ROI
% INPUTS:
%   roi [str] - name of the ROI
%       can be specific or general (e.g., DLPFC or LPFC)
% OUTPUTS:
%   RGB [3 floats] - RGB for this ROI
%
% Values taken from: https://www.w3schools.com/colors/colors_groups.asp

switch roi
    % General ROIs
    case 'LPFC'         % Orange
        RGB = [253 192 134]./256;        
    case 'MPFC'         % Purple
        RGB = [190 174 212]./256;
    case 'INS'          % Green
        RGB = [127 201 127]./256;
    case 'OFC'          % Dark Blue
        RGB = [56 108 176]./256;

    % LPFC Subregions - reds
    case 'FPC'
        RGB = [179 0 0]./256;
    case 'DLPFC'
        RGB = [227 74 51]./256;
    case 'VLPFC'
        RGB = [252 141 89]./256;
    case 'PM'
        RGB = [253 187 132]./256;
    case 'M1'
        RGB = [253 212 158]./256;
    case 'S1'
        RGB = [254 240 217]./256;
    
    % MPFC Subregions - blues
    case 'ACC'
        RGB = [8 81 156]./256;
    case 'preSMA'
        RGB = [49 130 189]./256;
    case 'aMCC'
        RGB = [107 174 214]./256;
    case 'SMA'
        RGB = [158 202 225]./256;
    case 'pMCC'
        RGB = [198 219 239]./256;
    case 'Precuneus'
        RGB = [239 243 255]./256;

    % Insula Subregions - greens
    case 'vaINS'
        RGB = [0 109 44]./256;
    case 'daINS'
        RGB = [49 163 84]./256;
    case 'FO'
        RGB = [116 196 118]./256;
    case 'mINS'
        RGB = [161 217 155]./256;
    case 'pINS'
        RGB = [199 233 192]./256;

    % OFC Subregions - yellows
    case 'mOFC'
        RGB = [255 215 0]./256;
    case 'lOFC'
        RGB = [255 255 0]./256;

    % Weird Cases
    case 'FWM'          % Gray
        RGB = [0.5 0.5 0.5];
    case ''
        RGB = [0.5 0.5 0.5];
    case 'OUT'          % Black
        disp('WARNING: Why are you trying to plot data mainly out of the brain???');
        RGB = [0 0 0];
end