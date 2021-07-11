function colors = fn_venn_colors(n_categories, varargin)
%% function colors = fn_venn_colors(n_categories)
% returns [n_cat, n_cat] matrix of colors for the number of categories
% Color determination:
%   if cond_id or model_lab is provided via varargin, pulls those
%       overlap is mean of cond_colors (or reg_colors)
%   if no cond_id:
%       n=3 colors based on R, G, B (Y, P, C)
%       n=2 colors based on G, B (C)
% Stat color matrix:
%   Primary: Red, Green, Blue
%   Combinations: Yellow, Pink, Cyan
%   R   Y   P
%   Y   G   C
%   P   C   B

%% Handle variables
% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'cond_id')
            cond_id = varargin{v+1};
        elseif strcmp(varargin{v},'model_lab')
            model_lab = varargin{v+1};
        elseif strcmp(varargin{v},'color_id')
            color_id = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

if ~exist('hard_colors','var'); color_id = 'def'; end

%% Determine Colors
if exist('cond_id','var')
    [~, cond_colors, ~] = fn_condition_label_styles(cond_id);
    colors = cell(numel(cond_colors));
    for c1 = 1:numel(cond_colors)
        for c2 = 1:numel(cond_colors)
            colors{c1,c2} = mean([cond_colors{c1};cond_colors{c2}],1);
        end
    end
elseif exist('model_lab','var')
    [reg_lab, ~, ~, reg_colors, ~] = fn_regressor_label_styles(model_lab);
    colors = cell(numel(reg_colors));
    if any(strcmp(model_lab,{'pWinPEus','ERPEs'}))
        colors{1,1} = [0 114 178]./256;     % EV - CB blue
        colors{2,2} = [204 121 167]./256;   % sRPE - CB reddish purple
        colors{3,3} = [0 158 115]./256;   % uRPE - CB bluish green
        colors{1,2} = [230 159 0]./256;      % EV + sRPE - CB orange
        colors{1,3} = [97 61 46]./256;   % EV + uRPE - brown
        colors{2,3} = [240 228 66]./256;   % sRPE + uRPE - CB yellow
        colors{3,1} = [0 0 0];              % all - black
        % EEG Communications Biology Colors:
%         colors{1,1} = [189 65 45]./256;     % EV - tomato (red)
%         colors{2,2} = [209 151 105]./256;   % sRPE - tan
%         colors{3,3} = [118 160 156]./256;   % uRPE - aqua
%         colors{1,2} = [97 61 46]./256;      % EV + sRPE - brown
%         colors{1,3} = [240 224 163]./256;   % EV + uRPE - creamy yellow
%         colors{2,3} = [188 147 196]./256;   % sRPE + uRPE - dark purple
%         colors{3,1} = [0 0 0];              % all - black
    elseif any(strcmp(model_lab,{'suRPE'}))
        colors{1,1} = [204 121 167]./256;   % sRPE - CB reddish purple
        colors{2,2} = [0 158 115]./256;     % uRPE - CB bluish green
        colors{1,2} = [240 228 66]./256;    % sRPE + uRPE - CB yellow
    else
        for r1 = 1:numel(reg_colors)
            for r2 = 1:numel(reg_colors)
                colors{r1,r2} = mean([reg_colors{r1};reg_colors{r2}],1);
            end
        end
    end
else
    % Define Colors
    if strcmp(color_id,'hard')
        r = {[1 0 0]./255};
        g = {[0 1 0]./255};
        b = {[0 0 1]./255};
        y = {[1 1 0]./255};
        p = {[1 0 1]./255};
        c = {[0 1 1]./255};
%     elseif strcmp('color_id','ROI')
%         r = {[255 158 163]./255};
%         g = {[163 255 158]./255};
%         b = {[158 163 255]./255};
%         y = {[255 250 146]./255};
%         p = {[250 146 255]./255};
%         c = {[146 255 250]./255};
    elseif strcmp(color_id,'def')
        r = {[255 158 163]./255};
        g = {[163 255 158]./255};
        b = {[158 163 255]./255};
        y = {[255 250 146]./255};
        p = {[250 146 255]./255};
        c = {[146 255 250]./255};
    else
        error(['unknown color_id: ' color_id]);
    end
    
    % Select Colors
    if n_categories==2
        colors = [g c; c b];
    elseif n_categories==3
        colors = [r y p; y g c; p c b];
    else
        error('n_categories must be 2 or 3');
    end
end

end