function condition_num = fn_condition_index(conditions, trl_info)
% Returns index of trial condition assignments based on requested conditions
% INPUTS:
%   conditions [str] - name of the set of conditions requested
%   trl_info [struct] - trial informatino structure contianing info for logic sorting
% OUTPUTS:
%   condition_n [int vector] - integer assignment of each trial based on conditions

[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);

condition_num = zeros(size(trl_info.trl_n));
for cond_ix = 1:numel(cond_lab)
    switch cond_lab{cond_ix}
        case 'Ez'
            condition_num(strcmp('easy',trl_info.cond)) = cond_ix;
        case 'Hd'
            condition_num(strcmp('hard',trl_info.cond)) = cond_ix;
        case 'Wn'
            condition_num(trl_info.hit==1) = cond_ix;
        case 'Ls'
            condition_num(trl_info.hit==0) = cond_ix;
        case 'Er'
            condition_num(trl_info.rt<trl_info.prdm.target) = cond_ix;
        case 'Lt'
            condition_num(trl_info.rt>trl_info.prdm.target) = cond_ix;
        case 'EzWn'
            matches = logical(strcmp('easy',trl_info.cond)) & logical(trl_info.hit==1);
            condition_num(matches) = cond_ix;
        case 'EzLs'
            matches = logical(strcmp('easy',trl_info.cond)) & logical(trl_info.hit==0);
            condition_num(matches) = cond_ix;
        case 'HdWn'
            matches = logical(strcmp('hard',trl_info.cond)) & logical(trl_info.hit==1);
            condition_num(matches) = cond_ix;
        case 'HdLs'
            matches = logical(strcmp('hard',trl_info.cond)) & logical(trl_info.hit==0);
            condition_num(matches) = cond_ix;
        otherwise
            error(['Invalid condition label: ' conditions]);
    end
end

if sum(condition_num==0)~=0
    warning(['Not all trials accounted for by conditions: ' conditions]);
end

end