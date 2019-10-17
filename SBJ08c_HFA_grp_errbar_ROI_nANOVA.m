function SBJ08c_HFA_grp_errbar_ROI_nANOVA(SBJs,stat_id,proc_id,an_id,atlas_id,roi_id,...
                                          plot_scat,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for active, RT correlation, and ANOVA epochs
%   RT correlation: any significance in stat_lim
%   ANOVA factors: any significance in stat_lim, after FDR correction
% OUTPUTS:
%   Bar chart with % ANOVA factors
% clear all; %close all;
% fig_ftype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end
if ischar(plot_scat); plot_scat = str2num(plot_scat); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir();
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Prep variables
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Get condition info
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

% Set up electrode counts
grp_cnt   = zeros([numel(SBJs) numel(roi_list) numel(grp_lab)]);
elec_cnt  = zeros([numel(SBJs) numel(roi_list)]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load ANOVA
    load([SBJ_vars.dirs.proc,SBJ,'_nANOVA_ROI_',stat_id,'_',an_id,'.mat']);
    
    %% Load ROI and GM/WM info
    load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    
%     % HACK!!!
%     cfgs = [];
%     if strcmp(SBJ,'CP24')
%         cfgs.channel = {'all','-RTO4'};
%     elseif strcmp(SBJ,'IR57')
%         cfgs.channel = {'all','-LAM4-5','-LAM5-6','-RIN8-9','-RIN9-10',...
%                         '-RSM1-2','-RSM2-3','-RSM3-4','-RTI9-10'};
%     elseif strcmp(SBJ,'IR68')
%         cfgs.channel = {'all','-LPC5-6','-LPC6-7','-LPC7-8'};
%     end
%     stat = ft_selectdata(cfgs,stat);

%     % Sort elecs by stat labels
%     if ~strcmp(SBJ,'IR66')  % HACK!!!
%         cfgs = []; cfgs.channel = stat.label;
%         elec = fn_select_elec(cfgs,elec);
%         elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
%     else
%         elec.roi = elec.(roi_id);
%     end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(w2.label)
        elec_ix = strcmp(w2.label{ch_ix},elec.label);
        if any(elec_ix)
            roi_ix = strcmp(elec.(roi_field){elec_ix},roi_list);
            if any(roi_ix)
                elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
                % Check for ANOVA group effects
                for grp_ix = 1:numel(grp_lab)
                    if any(w2.qval(grp_ix,ch_ix,:)<=0.05,3)
                        grp_cnt(sbj_ix,roi_ix,grp_ix) = grp_cnt(sbj_ix,roi_ix,grp_ix)+1;
                    end
                end
            else
                fprintf(2,'\t%s: %s in %s not in %s!\n',SBJ,w2.label{ch_ix},...
                    elec.(roi_field){elec_ix},roi_id);
            end
        else
            fprintf(2,'\t%s: %s in w2 missing from elec!\n',SBJ,w2.label{ch_ix});
        end
    end
    clear SBJ SBJ_vars w2 elec
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
if plot_scat
    scat_suffix = '_SBJscat';
else
    scat_suffix = '';
end
% Create and format the plot
fig_name = ['GRP_HFA_errbar_' stat_id '_' atlas_id '_' roi_id scat_suffix];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% Compile Data in Plotting Format
scat_vars = cell(size(grp_lab));
bar_data  = zeros([numel(grp_lab) numel(roi_list)]);
var_data  = zeros([numel(grp_lab) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
    for grp_ix = 1:numel(grp_lab)
        bar_data(grp_ix,roi_ix) = sum(grp_cnt(:,roi_ix,grp_ix))/sum(elec_cnt(:,roi_ix));
        var_data(grp_ix,roi_ix) = nanstd(grp_cnt(:,roi_ix,grp_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
        scat_vars{grp_ix} = squeeze(grp_cnt(:,:,grp_ix));
    end
end

% Plot Activations by ROI
b = cell(size(grp_lab));
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(roi_list));   %bar for activation, deactivation, condition
for grp_ix = 1:numel(grp_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{grp_ix} = bar(bar_offsets+grp_ix,diag(bar_data(grp_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(b{grp_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        % Plot individual subject percentages as scatters on top
        line([bar_offsets(roi_ix) bar_offsets(roi_ix)]+grp_ix,...
            [bar_data(grp_ix,roi_ix)+var_data(grp_ix,roi_ix) bar_data(grp_ix,roi_ix)-var_data(grp_ix,roi_ix)],...
            'Color','k','LineWidth',1.5);
        % Overlay scatter for individual SBJ
        if plot_scat
            has_elecs = squeeze(elec_cnt(:,roi_ix)~=0);
            s = scatter(scat_offsets(has_elecs)+bar_offsets(roi_ix)+grp_ix,...
                scat_vars{grp_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),50,'k*');
        end
    end
end
leg_loc = 'northwest';
if plot_scat
    legend([b{1},s],roi_list{:},'Individuals','Location',leg_loc);
else
    legend(b{1},roi_list{:},'Location',leg_loc);
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = [0.5 0.5+numel(grp_lab)];
ax.XTick   = 1:numel(grp_lab);
ax.XColor  = 'k';
ax.XTickLabel = grp_lab;

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLabel.FontSize = 16;
ax.YLim            = [0 1];%ax.YLim(2)];%[0 0.6];%
ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Proportion of Electrodes Showing Significant Effects';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
set(gca,'FontSize',16);

%% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP_errbar_ROI/'...
        stat_id '/' an_id '/' atlas_id '_' roi_id '/'];
    if ~exist(fig_dir,'dir')
        [~,~] = mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
