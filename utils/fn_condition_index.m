function condition_n = fn_condition_index(cond_lab, bhv)
% Returns index of trial condition assignments based on requested conditions
% INPUTS:
%   cond_lab [cell array] - strings with names of the set of conditions requested
%   bhv [struct] - trial information structure contianing info for logic sorting
%       fields: Total_Trial, Block, Feedback, RT, Timestamp, Tolerance, Trial, Hit, Score, bad_fb, Condition, ITI, ITI type
% OUTPUTS:
%   condition_n [int vector] - integer assignment of each trial based on conditions

condition_n = zeros(size(bhv.trl_n));
for cond_ix = 1:numel(cond_lab)
    switch cond_lab{cond_ix}
        % Target Time block and trial types
        case 'All'
            condition_n = ones(size(condition_n));
        case 'Ez'
            condition_n(strcmp('easy',bhv.cond)) = cond_ix;
        case 'Hd'
            condition_n(strcmp('hard',bhv.cond)) = cond_ix;
        case 'Wn'
            condition_n(strcmp(bhv.fb,'W')) = cond_ix;
        case 'Ls'
            condition_n(strcmp(bhv.fb,'L')) = cond_ix;
        case 'Nu'
            condition_n(strcmp(bhv.fb,'N')) = cond_ix;
        
        % Target Time specific conditions
        case 'EzWn'
            matches = logical(strcmp('easy',bhv.cond)) & strcmp(bhv.fb,'W');
            condition_n(matches) = cond_ix;
        case 'EzLs'
            matches = logical(strcmp('easy',bhv.cond)) & strcmp(bhv.fb,'L');
            condition_n(matches) = cond_ix;
        case 'EzNu'
            matches = logical(strcmp('easy',bhv.cond)) & strcmp(bhv.fb,'N');
            condition_n(matches) = cond_ix;
        case 'HdWn'
            matches = logical(strcmp('hard',bhv.cond)) & strcmp(bhv.fb,'W');
            condition_n(matches) = cond_ix;
        case 'HdLs'
            matches = logical(strcmp('hard',bhv.cond)) & strcmp(bhv.fb,'L');
            condition_n(matches) = cond_ix;
        case 'HdNu'
            matches = logical(strcmp('hard',bhv.cond)) & strcmp(bhv.fb,'N');
            condition_n(matches) = cond_ix;
%         case 'Ex'
%             match1 = logical(strcmp('easy',bhv.cond)) & logical(bhv.hit==1);
%             match2 = logical(strcmp('hard',bhv.cond)) & logical(bhv.hit==0);
%             condition_n(match1 | match2) = cond_ix;
%         case 'Ue'
%             match1 = logical(strcmp('easy',bhv.cond)) & logical(bhv.hit==0);
%             match2 = logical(strcmp('hard',bhv.cond)) & logical(bhv.hit==1);
%             condition_n(match1 | match2) = cond_ix;
            
        % Target Time conditions based on RPE valence
        case 'AllPos'
            matches = fn_condition_index({'EzWn','HdWn','HdNu'}, bhv);
            condition_n(matches~=0) = cond_ix;
        case 'AllNeg'
            matches = fn_condition_index({'EzLs','HdLs','EzNu'}, bhv);
            condition_n(matches~=0) = cond_ix;
        
        % Target Time RT performance assignment
        case 'Er'
            condition_n(bhv.rt<bhv.prdm.target) = cond_ix;
        case 'Lt'
            condition_n(bhv.rt>bhv.prdm.target) = cond_ix;
        
        % Oddball task conditions
        case 'Odd'
            matches = logical(strcmp('odd',bhv.cond));
            condition_n(matches) = cond_ix;
        case 'Std'
            matches = logical(strcmp('std',bhv.cond));
            condition_n(matches) = cond_ix;
        case 'Tar'
            matches = logical(strcmp('tar',bhv.cond));
            condition_n(matches) = cond_ix;
        otherwise
            error(['Invalid condition label: ' cond_lab{cond_ix}]);
    end
end

if sum(condition_n==0)~=0
    warning(['Not all trials accounted for by conditions: ' strjoin(cond_lab,',')]);
end

end