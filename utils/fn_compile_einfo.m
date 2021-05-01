function fn_compile_einfo(SBJ,proc_id)
%% Compile ROI info from single electrodes into bipolar pairs
%   Goes from smallest to largest (ELEC1-ELEC2, ELEC2-ELEC3, etc.)
%   Pairs are drawn from preprocessed data labels
%
% new_labels = cell array of strings with combined labels
%   can be used with ft_preprocessing as cfg_rereference.montage.labelnew
% weights = [length(new_labels),length(labels)] matrix of weights to combine elecs
%   can be used with ft_preprocessing as cfg_rereference.montage.tra
% left_out_ch = labels of channels that aren't included (don't have contiguous pair)

% Set up paths
addpath('/home/knight/hoycw/PRJ_Error/scripts/');
addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load data
eval(['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Error/scripts/proc_vars/' proc_id '_vars.m']);

% % Original (single electrode) labels
% import  = load([SBJ_vars.dirs.import SBJ '_1000hz.mat']);
% raw_lab = import.data.label;

% Bipolar ROI Labels
preproc = load([SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat']);
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
data = ft_selectdata(cfgs,preproc.data);
new_lab =  data.label;

%% Load ROI and GM/WM info
fprintf('\tReading roi csv file\n');
roi_file = fopen([SBJ_vars.dirs.SBJ SBJ '_elec_ROIs.tsv'], 'r');
% roi.csv contents:
%   Electrode, orig_ROI, Primary tissue, neighboring tissue, Gross/General ROI (gROI),
%   Detailed ROI, Anatomy Notes, Data notes
roi_info = textscan(roi_file, '%s %s %s %s %s %s %s %s', 'HeaderLines', 1,...
    'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
fclose(roi_file);
% WARNING: It seems to miss the last field (Data Notes) for the last row
% roi_list = unique([roi_info{:,6}]);
% roi_gen_list = unique([roi_info{:,5}]);

% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain
einfo = cell(numel(new_lab),7);
for e_ix = 1:numel(new_lab)
    % Electrode Names and Indices
    einfo{e_ix,1} = new_lab{e_ix};
    
    % CAR Reference
    if isempty(strfind(new_lab{e_ix},'-'))
        orig_e_ix = find(strcmp(roi_info{1},new_lab{e_ix}));
        einfo{e_ix,2} = roi_info{5}{orig_e_ix};     % ROI
        einfo{e_ix,3} = roi_info{4}{orig_e_ix};     % gROI
        einfo{e_ix,4} = '';                         % ROI2
        einfo{e_ix,5} = roi_info{2}{orig_e_ix};     % tissue
        if ~strcmp(einfo{e_ix,5},'GM')
            error(['ERROR: CAR referenced electrode (assumed to be grid) is not in GM! Check out: '...
                SBJ '_' new_lab{e_ix}]);
        end
        einfo{e_ix,6} = 1;                          % gm_weight
        einfo{e_ix,7} = 0;                          % out
    % Bipolar reference
    else
        orig_lab1 = new_lab{e_ix}(1:strfind(new_lab{e_ix},'-')-1);
        orig_e_ix1 = find(strcmp(roi_info{1},orig_lab1));
        probe = orig_lab1(regexp(orig_lab1,'\D')); % did this weird 2 line thing because MATLAB crashes if it's 1
        orig_lab2 = strcat(probe, new_lab{e_ix}(strfind(new_lab{e_ix},'-')+1:end));
        orig_e_ix2 = find(strcmp(roi_info{1},orig_lab2));
        if any([isempty(orig_e_ix1),isempty(orig_e_ix2)])
            disp(e_ix);
            error(['ERROR: Label not found. lab1=' orig_lab1 ' and lab2=' orig_lab2]);
        end
        
        % ROI, gROI, tissue, etc.
        %   [new_roi, new_groi, new_roi2, new_tissue, gm_weight, out] = ...
        [einfo{e_ix,2}, einfo{e_ix,3}, einfo{e_ix,4}, einfo{e_ix,5}, einfo{e_ix,6}, einfo{e_ix,7}] = ...
            fn_combine_ROI_bipolar_logic(roi_info{6}{orig_e_ix1}, roi_info{5}{orig_e_ix1},...    % ROI, gROI
            roi_info{3}{orig_e_ix1}, roi_info{4}{orig_e_ix1},...    % primary, secondary tissue
            roi_info{6}{orig_e_ix2}, roi_info{5}{orig_e_ix2},...    % ROI, gROI
            roi_info{3}{orig_e_ix2}, roi_info{4}{orig_e_ix2});      % primary, secondary tissue
    end
end

%% Save data
output_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_einfo_',proc_id,'.mat');
fprintf('============== Saving %s ==============\n',output_filename);
save(output_filename, '-v7.3', 'einfo');

end
