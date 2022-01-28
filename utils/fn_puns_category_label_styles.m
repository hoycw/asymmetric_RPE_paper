function [labels, names, colors, line_styles, markers] = fn_puns_category_label_styles(cat_id)
%% Converts the name of a category ID into category labels, plotting colors/styles
%   Assumes model_id = 'EpnRPEs_DifFB' ~ EV + pRPE + nRPE
%       RL Features: expected value (logistic regression), positive RPE, negative RPE
%   categories: (based on coefficient signs for pRPE and nRPE RL regressors)
%       pRPE:  positive RPE (increasing activity to more positive RPE)
%       nRPE:  negative RPE (increasing activity to more negative RPE)
%       sRPE:  signed RPE (increasing activity to more positive RPE, decreasing activity to more negative RPEs)
%       uRPE:  unsigned RPE/salience (increasing activity to larger RPEs regardless of valence)
%       ipRPE: inverse/inhibitory positive RPE (decreasing activity to more positive RPE)
%       inRPE: inverse/inhibitory negative RPE (decreasing activity to more negative RPE)
%       isRPE: inverse/inhibitory signed RPE (decreasing activity to more positive RPE, increasing activity to more negative RPEs)
%       iuRPE: inverse/inhibitory unsigned RPE/salience (decreasing activity to larger RPEs regardless of valence)
% INPUTS:
%   cat_id [str] - format is ['cat_label' '_' 'cond_lab']
% OUTPUTS:
%   labels [cell array] - string short-hand labels of specific categories
%   names [cell array] - string longer, full labels of specific categories
%   colors [cell array] - [R G B] tuples (scaled to 0-1) per category
%   line_styles [cell array] - line plotting styles per category
%   markers [cell array] - scatter plot markers per category
% Color options: Taken from colorbrewer2.org, qualitative, 5-class Set1
%   purple: [152 78 163]
%   orange: [255 127 0]
%   brown:  [166 86 40]
%   yellow: [255 255 51]
%   pink:   [247 129 191]
%   red:    [228 26 28]
%   green:  [55 126 184]
%   blue:   [77 175 74]
%   gray:   [0.5 0.5 0.5]
% Found in Online Image "color_combinations" on Desktop
%   Antiquated Aqua:    [118 160 156]
%   Warm Sunglow (tan): [209 151 105]
%   Tarrazzo Brown:     [97 61 46]
%   Tomato Tango:       [189 65 45]

%% List of possible regressors and their colors
regressors  = {...
    'pRPE','nRPE','sRPE','uRPE',... % main categories with excitatory coding
    'ipRPE','inRPE','isRPE','iuRPE'...  % inverse/inhibitory coding categories
    };
regressor_names = {...
    'Positive RPE','Negative RPE','Signed RPE','Unsigned RPE',...
    'inhib Positive RPE','inhib Negative RPE','inhib Signed RPE','inhib Unsigned RPE'...
    };
regressor_colors = {...
    [228 26 28]./255, [55 126 184]./255,... % colorbrewer red, colorbrewer blue
    [209 151 105]./255, [0 158 115]./255,... % CB reddish purple, CB bluish green
    [228 26 28]./255, [55 126 184]./255,... % colorbrewer red, colorbrewer blue
    [209 151 105]./255, [0 158 115]./255 ... % CB reddish purple, CB bluish green
    };
regressor_markers = {...
    'x','+','d','o',...
    'x','+','d','o' ...
    };

%% Convert cat_id into set of conditions
switch cat_id
    case 'puns'
        labels = {'pRPE','nRPE','sRPE','uRPE'};
    case 'ipuns'
        labels = {'pRPE','nRPE','sRPE','uRPE','ipRPE','inRPE','isRPE','iuRPE'};
        
    case 'pn'
        labels = {'pRPE','nRPE'};
    case 'ipn'
        labels = {'pRPE','nRPE','ipRPE','inRPE'};
        
    case 'su'
        labels = {'sRPE','uRPE'};
    case 'isu'
        labels = {'isRPE','iuRPE'};
        
    otherwise
        error(strcat('Unknown cat_id: ',cat_id));
end

%% Assign colors and line styles
names       = cell(size(labels));
colors      = cell(size(labels));
line_styles = cell(size(labels));
markers     = cell(size(labels));
for cat_ix = 1:numel(labels)
    names{cat_ix}   = regressor_names{strcmp(labels{cat_ix},regressors)};
    colors{cat_ix}  = regressor_colors{strcmp(labels{cat_ix},regressors)};
    markers{cat_ix} = regressor_markers{strcmp(labels{cat_ix},regressors)};
    if ~contains(labels{cat_ix},'i')
        line_styles{cat_ix} = '-';
    else
        line_styles{cat_ix} = ':';
    end
end

end
