function elec = fn_combine_ROI_bipolar_logic(reref,report_fname,atlas)
% Logic to determine the properties assigned to a bipolar referenced pair of electrodes
%   If already manually adjusted inputs, skips most comparison logic and
%   just flags any non-matching neighbors for inspection.
% INPUTS: 
%   reref [struct] - output elec struct of ft_apply_montage
%       .chanpos & .label should have new dimensions
%       .hemi, .atlas_lab, .atlas_lab2, .tissue, .tissue_prob have old dimensions and need updates
%           .tissue2 also has old labels, but for manual adjusted will not be updated
%   report_fname [str] - filename to write results
%   atlas [struct] - output of fn_load_recon_atlas (ft_read_atlas)
%       used to recompute atlas labels for averaged bipolar positions in case of a tie
% OUTPUTS:
%   elec [struct] - modified elec to adapt above fields for bipolar logic
%   results.txt [log file] - prints the outcome for each bipolar pair,
%   along with warnings for roi_flags

%% Set up log
log_fname = [report_fname '_log.txt'];
log = fopen(log_fname,'a');
% fprintf('Opening log file %s\n',log_fname);
fprintf(log,'%s\n',datestr(datetime));

%% Create initial copy
elec.chanpos       = reref.chanpos;
elec.coordsys      = reref.coordsys;
elec.elecpos       = reref.elecpos;
elec.label         = reref.label;
elec.unit          = reref.unit;
elec.type          = reref.type;
elec.atlas_id      = reref.atlas_id;
elec.tissue_labels = reref.tissue_labels;
elec.tra           = reref.tra;
elec.chanposold    = reref.chanposold;
elec.labelold      = reref.labelold;
elec.reref         = 1;

%% Combine pair info
elec.hemi         = cell(size(elec.label));
elec.atlas_lab    = cell(size(elec.label));
elec.atlas_prob   = zeros(size(elec.label));
elec.atlas_qryrng = zeros(size(elec.label));
elec.atlas_lab2   = cell(size(elec.label));
elec.tissue       = cell(size(elec.label));
% elec.tissue2      = cell(size(elec.label));
elec.gm_weight    = zeros(size(elec.label));
elec.tissue_prob  = zeros([numel(elec.label) numel(elec.tissue_labels)]);
elec.roi_flag     = zeros(size(elec.label));             % to flag problematic cases
elec.inputs       = cell(size(elec.label));              % keep the bipolar pair info to assess algorithm
if isfield(reref,'man_adj')
    manual = 1;
    elec.gROI       = cell(size(elec.label));
    elec.ROI        = cell(size(elec.label));
    elec.man_adj    = zeros(size(elec.label));
    elec.par_vol    = zeros(size(elec.label));
else
    manual = 0;
end
for p = 1:numel(elec.label)
    pair_ix = find(elec.tra(p,:)~=0);
    
    %% Tissue Logic
    elec.inputs{p}.tissue_prob = reref.tissue_prob(pair_ix,:);
    elec.inputs{p}.tissue      = reref.tissue(pair_ix);
    
    elec.tissue_prob(p,:) = mean(reref.tissue_prob(pair_ix,:),1);
    if any(strcmp(reref.tissue(pair_ix),'GM'))
        elec.tissue{p} = 'GM';
    else
        if any(~strcmp(reref.tissue(pair_ix),'WM'))
            fprintf(2,'\tWARNING: bad tissues %s and %s for elec %s!\n',reref.tissue{pair_ix},elec.label{p});
        end
        elec.tissue{p} = 'WM';
    end
%     [~,tiss_idx] = sort(elec.tissue_prob(p,:),'descend');
%     elec.tissue(p)  = elec.tissue_labels(tiss_idx(1));
%     elec.tissue2(p) = elec.tissue_labels(tiss_idx(2));
    
    % Manual computation
    elec.inputs{p}.gm_weight = [reref.gm_weight(pair_ix(1)) reref.gm_weight(pair_ix(2))];
    elec.gm_weight(p) = mean(elec.inputs{p}.gm_weight);
    
    % Manual adjustment fields
    if manual
        elec.inputs{p}.ROI        = [reref.ROI(pair_ix(1)) reref.ROI(pair_ix(2))];
        elec.inputs{p}.gROI       = [reref.gROI(pair_ix(1)) reref.gROI(pair_ix(2))];
        elec.inputs{p}.man_adj    = [reref.man_adj(pair_ix(1)) reref.man_adj(pair_ix(2))];
        elec.inputs{p}.par_vol    = [reref.par_vol(pair_ix(1)) reref.par_vol(pair_ix(2))];
        elec.inputs{p}.anat_notes = [reref.anat_notes(pair_ix(1)) reref.anat_notes(pair_ix(2))];
    end
    
    %% ROI Logic
    % 1. GM primary label
    %   1a. match a secondary label?
    %   1b. greater GM tissue_prob?
    %   1c. tie = choose ROI at middle position between orig elecpos
    % 2. Greater GM tissue_prob
    elec.inputs{p}.atlas_lab    = [reref.atlas_lab(pair_ix(1)) reref.atlas_lab(pair_ix(2))];
    elec.inputs{p}.atlas_prob   = [reref.atlas_prob(pair_ix(1)) reref.atlas_prob(pair_ix(2))];
    elec.inputs{p}.atlas_qryrng = [reref.atlas_qryrng(pair_ix(1)) reref.atlas_qryrng(pair_ix(2))];
    elec.inputs{p}.atlas_lab2   = [reref.atlas_lab2(pair_ix(1)) reref.atlas_lab2(pair_ix(2))];
    elec.inputs{p}.atlas_prob2  = [reref.atlas_prob2(pair_ix(1)) reref.atlas_prob2(pair_ix(2))];
    
    % Select based on GM label
    non_gm = 0;
    gm_ix  = strcmp(reref.tissue_labels,'GM');
    if manual
        elec.roi_flag(p) = 1;
        if strcmp(reref.ROI{pair_ix(1)},reref.ROI{pair_ix(2)})  % Same ROI
            roi_ix = 1;
            elec.roi_flag(p) = 0;
        elseif all(strcmp({reref.tissue{pair_ix(1)}, reref.tissue{pair_ix(2)}},'GM'))   % Both GM
            % Check if only one is partial volume (pick other one)
            if reref.par_vol(pair_ix(1)) && ~all(reref.par_vol(pair_ix))
                roi_ix = 2;
            elseif reref.par_vol(pair_ix(2)) && ~all(reref.par_vol(pair_ix))
                roi_ix = 1;
            % Compare gm_weights
            elseif reref.gm_weight(pair_ix(1))~=reref.gm_weight(pair_ix(2))
                [~,roi_ix] = max(reref.gm_weight(pair_ix));
            else
                roi_ix = 1; % Arbitrarily pick most medial, then fix later!
                fprintf(2,'WARNING: GM Tie for %s between %s and %s!\n',elec.label{p},...
                                    reref.ROI{pair_ix(1)}, reref.ROI{pair_ix(2)});
                fprintf(log,'WARNING: GM Tie for %s between %s and %s!\n',elec.label{p},...
                                    reref.ROI{pair_ix(1)}, reref.ROI{pair_ix(2)});
            end
        elseif strcmp(reref.tissue{pair_ix(1)},'GM')    % One is GM
            roi_ix = 1;
        elseif strcmp(reref.tissue{pair_ix(2)},'GM')    % One is GM
            roi_ix = 2;
        else                                            % Neither GM
%             non_gm = 1;
            % Check if one has WM as gROI/ROI
            if ~strcmp(reref.gROI{pair_ix(1)},'WM') && strcmp(reref.gROI{pair_ix(2)},'WM')
                roi_ix = 1;
            elseif strcmp(reref.gROI{pair_ix(1)},'WM') && ~strcmp(reref.gROI{pair_ix(2)},'WM')
                roi_ix = 2;
            % Check if only one is on edge of GM
            elseif ~isempty(strfind(reref.anat_notes{pair_ix(1)},'edge')) && isempty(strfind(reref.anat_notes{pair_ix(2)},'edge'))
                roi_ix = 1;
            elseif isempty(strfind(reref.anat_notes{pair_ix(1)},'edge')) && ~isempty(strfind(reref.anat_notes{pair_ix(2)},'edge'))
                roi_ix = 2;
            else
                % Check if only one is partial volume (pick that one)
                if reref.par_vol(pair_ix(1)) && ~reref.par_vol(pair_ix(2))
                    roi_ix = 1;
                elseif ~reref.par_vol(pair_ix(1)) && reref.par_vol(pair_ix(2))
                    roi_ix = 2;
                % Compare gm_weights
                elseif reref.gm_weight(pair_ix(1))~=reref.gm_weight(pair_ix(2))
                    [~,roi_ix] = max(reref.gm_weight(pair_ix));
                % Competitive logic: Greater GM?
                elseif reref.tissue_prob(pair_ix(1),gm_ix) > reref.tissue_prob(pair_ix(2),gm_ix)
                    roi_ix = 1;
                elseif reref.tissue_prob(pair_ix(1),gm_ix) < reref.tissue_prob(pair_ix(2),gm_ix)
                    roi_ix = 2;
                else
                    roi_ix = 1; % Arbitrarily pick most medial, then fix later!
                    fprintf(2,'WARNING: non-GM Tie for % s between %s and %s!\n',elec.label{p},...
                        reref.ROI{pair_ix(1)}, reref.ROI{pair_ix(2)});
                    fprintf(log,'WARNING: non-GM Tie for % s between %s and %s!\n',elec.label{p},...
                        reref.ROI{pair_ix(1)}, reref.ROI{pair_ix(2)});
                end
            end
        end
    elseif all(strcmp({reref.tissue{pair_ix(1)}, reref.tissue{pair_ix(2)}},'GM'))   % Both GM
        % Start using the ROIs identified using the optimal query range
        if strcmp(reref.atlas_lab{pair_ix(1)},reref.atlas_lab{pair_ix(2)}) || any(strcmp(reref.atlas_lab{pair_ix(1)},reref.atlas_lab2{pair_ix(2)}))
            % Both are the same ROI, or lab(1) matches secondary lab2(2)
            roi_ix = 1;
        elseif any(strcmp(reref.atlas_lab2{pair_ix(1)},reref.atlas_lab{pair_ix(2)}))
            % lab(2) matches a secondary of lab2(1)
            roi_ix = 2;
        else
            % Competitive logic: Greater GM?
            if reref.tissue_prob(pair_ix(1),gm_ix) < reref.tissue_prob(pair_ix(2),gm_ix)
                roi_ix = 2;
            elseif  reref.tissue_prob(pair_ix(1),gm_ix) > reref.tissue_prob(pair_ix(2),gm_ix)
                roi_ix = 1;
            else    % TIE! choose whichever matches the middle position
                fprintf(2,'WARNING: GM Tie for %s between %s and %s!\n',elec.label{p},...
                                    reref.atlas_lab{pair_ix(1)}, reref.atlas_lab{pair_ix(2)});
                roi_ix = 3;
                elec.roi_flag(p) = 1;
            end
        end
    elseif strcmp(reref.tissue{pair_ix(1)},'GM')    % One is GM
        roi_ix = 1;
    elseif strcmp(reref.tissue{pair_ix(2)},'GM')    % One is GM
        roi_ix = 2;
    else                                            % Neither GM
        non_gm = 1;
        % Competitive logic: Greater GM?
        if reref.tissue_prob(pair_ix(1),gm_ix) < reref.tissue_prob(pair_ix(2),gm_ix)
            roi_ix = 2;
        elseif reref.tissue_prob(pair_ix(1),gm_ix) > reref.tissue_prob(pair_ix(2),gm_ix)
            roi_ix = 1;
        else    % TIE! choose whichever matches the middle position
            fprintf(2,'WARNING: Non-GM TIE for %s between %s and %s!\n',elec.label{p},...
                reref.atlas_lab{pair_ix(1)}, reref.atlas_lab{pair_ix(2)});
            roi_ix = 3;
            elec.roi_flag(p) = 1;
        end
    end
    
    %% Check hemisphere match
    if ~strcmp(reref.hemi{pair_ix(1)},reref.hemi{pair_ix(2)})
        warning(['Hemispheres do not match between ' reref.label{pair_ix(1)} ' and ' reref.label{pair_ix(2)}]);
        elec.roi_flag(p) = 1;
    end
    
    %% Add selected ROIs or atlas_labels
    if roi_ix==3
        % TIE
        if ~strcmp(reref.hemi{pair_ix(1)},reref.hemi{pair_ix(2)})
            error('Should not be hemisphere mismatch between GM_prob tie!');
        end
        elec.hemi{p} = reref.hemi{pair_ix(1)};
        
        % For tie, choose the averaged position ROI
        cfgs.channel = elec.label(p);
        if numel(reref.label) > 1
            elec_tmp = fn_select_elec(cfgs,reref);
        else
            elec_tmp = reref;
        end
        elec_tmp = fn_atlas_lookup(elec_tmp,atlas,'min_qry_rng',1,'max_qry_rng',5);
        elec.atlas_lab{p}  = elec_tmp.atlas_lab{1};
        elec.atlas_prob(p) = elec_tmp.atlas_prob(1);
        elec.atlas_qryrng(p) = elec_tmp.atlas_qryrng(1);
        if ~isempty(elec_tmp.atlas_lab2{1})
            elec.atlas_lab2{p}  = elec_tmp.atlas_lab2{1}{1};
        else
            elec.atlas_lab2{p}  = '';
        end
        
        fprintf(log,'TIE WARNING: %s (%s vs. %s) = %s\n',elec.label{p},...
            reref.atlas_lab{pair_ix(1)}, reref.atlas_lab{pair_ix(2)}, elec.atlas_lab{p});
    elseif non_gm
        % Both non-GM
        elec.hemi{p}       = reref.hemi{pair_ix(roi_ix)};
        % Search for GM labels
        new_labs = cell([2 1]); new_probs = zeros([2 1]);
        for i = 1:2
            % Select 1st/2nd labels in order of which elec won logic above
            if i==1; sel_ix = roi_ix; else sel_ix = roi_ix~=[1 2]; end
            
            % Check for secondary labels that are GM
            gm_match = strcmp(fn_atlas2roi_labels(reref.atlas_lab2{pair_ix(sel_ix)},reref.atlas_id,'tissue'),'GM');
            if any(gm_match)
                % If GM secondary label found, choose most likely GM label
                [~, prob_idx] = sort(reref.atlas_prob2{pair_ix(sel_ix)},'descend');
                % Narrow to only GM secondary labels
                prob_idx = prob_idx(gm_match);
                new_labs{i}  = reref.atlas_lab2{pair_ix(sel_ix)}{prob_idx(1)};
                new_probs(i) = reref.atlas_prob2{pair_ix(sel_ix)}(prob_idx(1));
            else
                % If no GM secondary label found, go with main label
                new_labs{i}  = reref.atlas_lab{pair_ix(sel_ix)};
                new_probs(i) = reref.atlas_prob(pair_ix(sel_ix));
            end
        end
        elec.atlas_lab{p}    = new_labs{1};
        elec.atlas_prob(p)   = new_probs(1);
        elec.atlas_qryrng(p) = reref.atlas_qryrng(pair_ix(roi_ix));
        elec.atlas_lab2{p}   = new_labs{2};
    else
        % GM Label
        elec.hemi{p}         = reref.hemi{pair_ix(roi_ix)};
        elec.atlas_lab{p}    = reref.atlas_lab{pair_ix(roi_ix)};
        elec.atlas_prob(p)   = reref.atlas_prob(pair_ix(roi_ix));
        elec.atlas_qryrng(p) = reref.atlas_qryrng(pair_ix(roi_ix));
        elec.atlas_lab2{p}   = reref.atlas_lab{pair_ix(roi_ix~=[1 2])};
        
        % Manually adjusted fields
        if manual
            elec.ROI{p}     = reref.ROI{pair_ix(roi_ix)};
            elec.gROI{p}    = reref.gROI{pair_ix(roi_ix)};
            elec.man_adj(p) = mean([reref.man_adj(pair_ix(1)) reref.man_adj(pair_ix(2))]);
            elec.par_vol(p) = mean([reref.par_vol(pair_ix(1)) reref.par_vol(pair_ix(2))]);
            % Concatenate notes while documenting which came for which
            note_str = cell([2 1]);
            for e = 1:2
                note_str{e} = elec.inputs{p}.anat_notes{e};
                if ~isempty(elec.inputs{p}.anat_notes{e})
                    note_str{e} = [num2str(e) ': ' note_str{e}];
                end
            end
            % Add separter if both have str
            if ~isempty(note_str{1}) && ~isempty(note_str{2})
                note_sep = '; ';
            else
                note_sep = '';
            end
            elec.anat_notes{p} = [note_str{1} note_sep note_str{2}];
        end
    end
    
    %% Log Results
    if manual
        fprintf(log,'%s = %s (%s | %s) [%.2f %.2f %.2f %.2f] vs %.2f\n',...
            elec.label{p}, elec.ROI{p}, elec.inputs{p}.ROI{1}, elec.inputs{p}.ROI{2},...
            elec.tissue_prob(p,1), elec.tissue_prob(p,2), elec.tissue_prob(p,3), elec.tissue_prob(p,4),...
            elec.gm_weight(p));
    else
        fprintf(log,'%s = %s (%.3f | %s) [%.2f %.2f %.2f %.2f] vs %.2f; %s (%.2f) and %s (%.2f)\n',...
            elec.label{p}, elec.atlas_lab{p}, elec.atlas_prob(p), elec.atlas_lab2{p},...
            elec.tissue_prob(p,1), elec.tissue_prob(p,2), elec.tissue_prob(p,3), elec.tissue_prob(p,4),...
            elec.gm_weight(p), elec.inputs{p}.atlas_lab{1}, elec.inputs{p}.atlas_prob(1),...
            elec.inputs{p}.atlas_lab{2}, elec.inputs{p}.atlas_prob(2));
    end
end

fclose(log);

end