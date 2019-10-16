function roi_lab = fn_atlas_roi_select_mesh(atlas_id, roi_id, hemi)
%% Returns ROI labels for plotting that ROI's mesh in a given atlas

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
tsv_fname = [root_dir 'PRJ_Error/data/atlases/atlas_mappings/atlas_ROI_mappings_' atlas_id '_both_hemi.tsv'];
fprintf('\tReading roi csv file: %s\n', tsv_fname);
roi_file = fopen(tsv_fname, 'r');
% roi.csv contents:
%   atlas_label, gROI_label, ROI_label, tissue, Known Variability, Alternative ROIs, Notes (not used)
roi_map = textscan(roi_file, '%s %s %s %s %d %s %s', 'HeaderLines', 1,...
    'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
fclose(roi_file);


roi_lab = {};
switch roi_id
    case 'lat'
        roi_lab = [...
            roi_map{1}(strcmp(roi_map{2},'LPFC'));
            roi_map{1}(strcmp(roi_map{2},'SM'));
            roi_map{1}(strcmp(roi_map{2},'OFC'));
            roi_map{1}(strcmp(roi_map{2},'TMP'));
            roi_map{1}(strcmp(roi_map{2},'OCC'));
            roi_map{1}(strcmp(roi_map{3},'SPL'));
            roi_map{1}(strcmp(roi_map{3},'IPL'))];
    case 'deep'
        roi_lab = [...
            roi_map{1}(strcmp(roi_map{2},'INS'))];
            %roi_map{1}(strcmp(roi_map{3},'AMG'));
            %roi_map{1}(strcmp(roi_map{3},'HPC'))];
    case 'LPFC'
        roi_lab = [...
            roi_map{1}(strcmp(roi_map{2},'LPFC'));
            roi_map{1}(strcmp(roi_map{2},'OFC'));
            'G_postcentral'; 'S_postcentral';           % PAR
            'G_front_sup'; 'G_and_S_paracentral'];      % MPFC
%             roi_map{1}(strcmp(roi_map{2},'INS'));
    case 'MPFC'
        roi_lab = [...
            roi_map{1}(strcmp(roi_map{2},'MPFC'));
            'G_cingul-Post-dorsal';                                     % PAR
            'S_cingul-Marginalis'; 'G_and_S_paracentral';               % SM
            'G_subcallosal'; 'G_rectus'; 'S_suborbital';                % OFC
            'G_front_sup'; 'G_and_S_transv_frontopol'];                 % LPFC
%             'S_subparietal'; 'G_precuneus'; 'G_cingul-Post-ventral';    % PAR
    case 'INS'
        roi_lab = roi_map{1}(strcmp(roi_map{2},'INS'));
    case 'OFC'
        roi_lab = roi_map{1}(strcmp(roi_map{2},'OFC'));
    case 'PAR'
        roi_lab = roi_map{1}(strcmp(roi_map{2},'PAR'));
    case 'TMP'
        roi_lab = roi_map{1}(strcmp(roi_map{2},'TMP'));
end

% Add hemisphere
if strcmp(hemi,'b')
    roi_lab = [roi_lab roi_lab];
end
for l = 1:size(roi_lab,1)
    if strcmp(hemi,'b')
        if any(strcmp(roi_lab(l,:),'Amygdala')) || any(strcmp(roi_lab(l,:),'Hippocampus'))
            if strcmp(hemi,'l')
                roi_lab{l,1} = ['Left-' roi_lab{l,1}];
            else
                roi_lab{l,2} = ['Right-' roi_lab{l,2}];
            end
        else
            roi_lab{l,1} = ['ctx_lh_' roi_lab{l,1}];
            roi_lab{l,2} = ['ctx_rh_' roi_lab{l,2}];
        end
    else
        if any(strcmp(roi_lab(l,:),'Amygdala')) || any(strcmp(roi_lab(l,:),'Hippocampus'))
            if strcmp(hemi,'l')
                roi_lab{l} = ['Left-' roi_lab{l}];
            else
                roi_lab{l} = ['Right-' roi_lab{l}];
            end
        else
            roi_lab{l} = ['ctx_' hemi 'h_' roi_lab{l}];
        end
    end
end
if strcmp(hemi,'b')
    roi_lab = [roi_lab(:,1); roi_lab(:,2)];
end
% roi_lab = roi_lab(~cellfun(@isempty,strfind(roi_lab,['_' hemi 'h_'])));
% atlas_labels = roi_map{1}(strcmp(roi_map{2},roi_id));
% % Select hemisphere
% atlas_labels = atlas_labels(~cellfun(@isempty,strfind(atlas_labels,['_' hemi 'h_'])));


end
