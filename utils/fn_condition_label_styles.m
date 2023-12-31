function [labels, names, colors, line_styles, markers] = fn_condition_label_styles(grp_id)
%% Converts the name of a set of conditions into labels, plotting colors/styles
%   Mostly for Target Time, but also some Oddball conditions
% INPUTS:
%   grp_id [str] - label for group of conditions to select
%       'Dif': {'Ez','Hd'} (difficulty)
%       'FB': {'Wn','Nu','Ls'} (feedback outcomes)
%       'DifFB': {'EzWn','EzNu','EzLs','HdWn','HdNu','HdLs'} (all outcomes)
%           NOTE: this should be default for selecting all trials/conditions;
%           'All' should only be used when intending to average across all conditions
% OUTPUTS:
%   labels [cell array] - string short-hand labels of specific conditions
%   names [cell array] - string longer, full labels of specific conditions
%   colors [cell array] - [R G B] tuples (scaled to 0-1) per condition
%   line_styles [cell array] - solid ('-') or dotted ('--') per condition
%   markers [cell array] - {'*','o','d'} markers per condition
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%   light red: [251 154 153]
%   dark red: [227 26 28]
%   light blue: [166 206 227]
%   dark blue: [31 120 180]
%   light green: [178 223 138]
%   dark green: [51 160 44]

%% List of possible labels and their colors
% Condition codes
conditions  = {...
    'Std','Tar','Odd',...                           % Oddball conditions
    'Ez','Hd','Wn','Ls','Nu',...                    % Target Time block and outcome types
    'EzWn','EzLs','EzNu','HdWn','HdLs','HdNu',...   % Target Time conditions
    'All','AllPos','AllNeg',...                     % Combinations of Target Time conditions
    'Er','Lt','ErQ1','ErQ2','MdQ3','LtQ4','LtQ5'};  % Performance-based target time selections

% Condition labels (longer strings for clear plotting legends)
condition_names = {...
    'Standard','Target','Novel',...
    'Easy','Hard','Win','Loss','Neutral',...
    'Easy Win','Easy Loss','Easy Neutral','Hard Win','Hard Loss','Hard Neutral',...
    'All Trials', 'All Positive', 'All Negative',...
    'Early','Late','Early Q1','Early Q2','Middle Q3','Late Q4','Late Q5'};

% Condition Colors ([R G B] scaled to 0-1)
condition_colors = {...
    [168 180 165]./256, [56 108 176]./256, [197,27,138]./256, ...
    [228,26,28]./256, [55,126,184]./256, [51 160 44]./256, [227 26 28]./256, [31 120 180]./256, ...
    [51 160 44]./256, [227 26 28]./256, [31 120 180]./256, ... % darker colors for easy 
    [51 160 44]./256, [227 26 28]./256, [31 120 180]./256, ... % darker colors for hard too 
    [0 0 0], [31 120 180]./256, [227 26 28]./256,...
    [166 97 26]./256, [1 133 113]./256, ... % brown, aqua
    [166 97 26]./256, [223 194 125]./256, [0.4 0.4 0.4], [128 205 193]./256, [1 133 113]./256 ...    % [brown, tan, gray, teal, aqua]
    };
% lighter RGB for loss/win/neutral:
%     [178 223 138]./256, [251 154 153]./256, [166 206 227]./256, ... % lighter colors for easy
% Colorblind friendly:
%     [0, 158, 115]./256, [213, 94, 0]./256, [0, 114, 178]./256, ... % easy: Vermillion, blue, bluish green
%     [240, 228, 66]./256, [204, 121, 167]./256, [86, 180, 233]./256, ... % hard: reddish purple, sky blue, yellow
% OLD Oddballs:
%   Target = [144 205 229]./256 light blue
%   Novel/Odd = [142 82 126]./256 purple

%% Convert grp_id into set of conditions
switch grp_id
    % ---------------------------------------------------------------------
    % Oddball Conditions
    case 'OB'
        labels = {'Std','Tar','Odd'};
    case 'rare'
        labels = {'Odd', 'Tar'};
    case 'Odd'
        labels = {'Odd'};
    case 'Tar'
        labels = {'Tar'};
    
    % ---------------------------------------------------------------------
    % Combinations of Target Time Conditions
    case 'All'
        labels = {'All'};
    case 'AllNeg'
        labels = {'AllNeg'};
    case 'AllPos'
        labels = {'AllPos'};
    case 'AllLrg'
        labels = {'AllLrg'};
    case 'Dif'                          % Difficulty
        labels = {'Ez', 'Hd'};
    case {'OutS','FB'}                  % Feedback (includes neutral)
        labels = {'Wn', 'Ls', 'Nu'};
    case 'Out'                          % Outcome (no neutral)
        labels = {'Wn', 'Ls'};
    case {'Val','RewP'}                          % Valence
        labels = {'AllPos', 'AllNeg'};
    
    % Collections of individual Target Time conditions
    case {'DifFB','Pos-Neg'}          % All 6 main conditions, no longer 'DifOutSur'
        labels = {'EzWn', 'EzNu', 'EzLs', 'HdWn', 'HdNu', 'HdLs'};
    case 'EzNu'                         % Neutral Outcomes in Easy
        labels = {'EzNu'};
    case 'HdNu'                         % Neutral Outcomes in Hard
        labels = {'HdNu'};
    case 'EHNu'                         % Neutral Outcomes in Easy+Hard
        labels = {'EzNu','HdNu'};
    case {'DifOut','DifOutUE','DifOutWL','DifOutdO','DifOutDO','Holroyd'}   % Four main conditions
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs'};
    case {'DifOutS', 'DifOutSx'}        % Four main conditions + all neutral outcomes combined
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs', 'Nu'};
%     case 'Dif*Out'
%         labels = {'Ex', 'Ue'};
%         colors = {[1,0,1]./256, [0,1,1]./256};  % magenta and cyan for now
    case 'EzOut'                        % Win/Loss in Easy
        labels = {'EzWn', 'EzLs'};
    case 'EzOutS'                       % All outcomes in Easy
        labels = {'EzWn', 'EzLs', 'EzNu'};
    case 'HdOut'                        % Win/Loss in Hard
        labels = {'HdWn', 'HdLs'};
    case 'HdOutS'                       % All outcomes in Hard
        labels = {'HdWn', 'HdLs', 'HdNu'};
    case 'Pos'                          % Conditions with positive RPE valence
        labels = {'EzWn','HdWn','HdNu'};
    case 'Neg'                          % Conditions with negative RPE valence
        labels = {'EzLs','HdLs','EzNu'};
    case 'Large-Small'                  % Conditions matched for valence and probability
        labels = {'EzNu','EzLs','HdWn','HdNu'};
    case 'Unlik-Lik'                    % Conditions matched for valence and magnitude
        labels = {'EzWn','EzNu','HdNu','HdLs'};
    
    % Performance (RT) based trial selection
    case 'Tim'                          % Early/Late RTs
        labels = {'Er', 'Lt'};
%     case 'Tar2'                         % Performance split of early/late
%         labels = {'Er', 'Lt'};
%     case 'Tar5'                         % RT-based split into quintiles
%         labels = {'ErQ1','ErQ2','MdQ3','LtQ4','LtQ5'};
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',grp_id));
end

%% Assign colors and line styles
names       = cell(size(labels));
colors      = cell(size(labels));
line_styles = cell(size(labels));
markers     = cell(size(labels));
for cond_ix = 1:numel(labels)
    names{cond_ix}  = condition_names{strcmp(labels{cond_ix},conditions)};
    colors{cond_ix} = condition_colors{strcmp(labels{cond_ix},conditions)};
    
    % Define Line Styles
    if isempty(strfind(labels{cond_ix},'Hd'))
        line_styles{cond_ix} = '-';     % Easy
    else
        line_styles{cond_ix} = '--';%'-.';    % Hard
    end
    
    % Define Marker Styles
    if strcmp(labels{cond_ix},'Nu')
        markers{cond_ix} = '*';     % Neutral (difficulty agnostic)
    elseif isempty(strfind(labels{cond_ix},'Hd'))
        markers{cond_ix} = 'o';     % Easy
    else
        markers{cond_ix} = 'd';    % Hard
    end
end

end
