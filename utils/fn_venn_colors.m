function colors = fn_venn_colors(n_categories, varargin)
%% function colors = fn_venn_colors(n_categories)
% returns [n_cat, n_cat] matrix of colors for the number of categories
% Color determination:
%   if cond_id is provided via varargin, pulls those
%       overlap is mean of cond_colors
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
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Determine Colors
if exist('cond_id','var')
    [~, cond_colors, ~] = fn_condition_label_styles(cond_id);
    colors = cell(numel(cond_colors));
    for c1 = 1:numel(cond_colors)
        for c2 = 1:numel(cond_colors)
            colors{c1,c2} = mean([cond_colors{c1};cond_colors{c2}],1);
        end
    end
else
    % Define Colors
    r = {[255 158 163]./255};
    g = {[163 255 158]./255};
    b = {[158 163 255]./255};
    y = {[255 250 146]./255};
    p = {[250 146 255]./255};
    c = {[146 255 250]./255};
    
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