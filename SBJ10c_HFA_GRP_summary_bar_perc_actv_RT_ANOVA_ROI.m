function SBJ10c_HFA_GRP_summary_bar_perc_actv_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,an_id,actv_win,roi_id,...
                                                            plt_id,save_fig,fig_vis,fig_filetype)
% Load HFA analysis results for active, RT correlation, and ANOVA epochs
%   Active- must be significant, and for S-locked, before mean(RT)
%   RT correlation: any significance in stat_lim
%   ANOVA factors: any significance in stat_lim, after FDR correction
% OUTPUTS:
%   Bar chart with % active, % deactivated, % RT correlations, % ANOVA factors
% clear all; %close all;
% fig_filetype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Error/scripts/');
addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Prep variables
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
% if rt_correlation
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
% end

% Load all ROI info
load('~/PRJ_Error/data/full_roi_lists.mat');
[roi_list, roi_colors, einfo_roi_col] = fn_roi_label_styles(roi_id);

% Get event timing
if strcmp(an_id(1:5),'HGm_S')
    event_lab = 'stim';
elseif strcmp(an_id(1:5),'HGm_R')
    event_lab = 'resp';
end

% Set up electrode counts
actv_cnt  = zeros([numel(SBJs) numel(roi_list)]);
dact_cnt  = zeros([numel(SBJs) numel(roi_list)]);
rt_cnt    = zeros([numel(SBJs) numel(roi_list)]);
grp_cnt   = zeros([numel(SBJs) numel(roi_list) numel(grp_lab)]);
elec_cnt  = zeros([numel(SBJs) numel(roi_list)]);

mean_RTs = NaN(size(SBJs));

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event_lab,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');
        mean_RTs(sbj_ix) = mean(trl_info.rt); % converts sec to ms
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
    load(actv_filename,'actv_ch','actv_ch_epochs');
    tmp = load(actv_filename,'hfa'); hfa_actv = tmp.hfa;
    
    %% Load ROI and GM/WM info
    einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
    load(einfo_filename);
    % Electrode Info Table:
    %   label- name of electrode
    %   ROI- specific region
    %   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
    %   ROI2- specific region of second electrode
    %   tissue- primary tissue type
    %   GM weight- percentage of electrode pair in GM
    %   Out- 0/1 flag for whether this is partially out of the brain
    
    % Sort by gROI, then ROI
    einfo = sortrows(einfo,[3,2]);
    if ~isempty(setdiff(stat.label,einfo(:,1)))
        error('ERROR: Electrodes do not match between stat and einfo!');
    end
    
    %% Process parameters
    %!!! is this the best way to do this??? Maybe not...
    sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
    
    % Restrict hfa_actv to stat_lim (should obviate some of below comparisons)
    cfg_lim = [];
    cfg_lim.latency = stat_lim;
    hfa_actv = ft_selectdata(cfg_lim,hfa_actv);
    hfa      = ft_selectdata(cfg_lim,hfa);
    stat     = ft_selectdata(cfg_lim,stat);
    
    % FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    % Confirm channels and time axis are the same
    %   All the rounding and such is because some stupid rounding errors...
    same_start = round(hfa_actv.time(1)*uint8(sample_rate))==round(stat.time(1)*uint8(sample_rate));
    same_end   = round(hfa_actv.time(end)*uint8(sample_rate))==round(stat.time(end)*uint8(sample_rate));
    same_numel = size(hfa_actv.time,2)==size(stat.time,2);
    if ~same_start || ~same_end || ~same_numel
        error('time axes are not the same across hfa analyses!');
    end
    if ~isempty(setdiff(stat.label,hfa_actv.label));
        error('Different electrodes across hfa analyses!');
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        both_ix = 0;
        einfo_ix = strmatch(stat.label(ch_ix),einfo(:,1),'exact');
        
        % If elec matches roi_list, get stats
        if ~isempty(strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact'))
            roi_ix = strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact');
            elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
            
            % Check if active, get epochs
            if ~isempty(strmatch(stat.label{ch_ix},actv_ch,'exact'))
                % Find significant epoch indices
                actv_ix = strmatch(stat.label{ch_ix},actv_ch,'exact');
                actv_epochs = actv_ch_epochs{actv_ix};
                
                % Toss late epochs if S-locked
                if strcmp('stim',event_lab)
                    actv_epochs(actv_epochs(:,1)<mean_RTs(sbj_ix),:) = [];
                end
                
                % Find sign of (de)activation
                actv_ep_sign = NaN([1 size(actv_epochs,1)]);
                sig_chunk_ix = NaN([1 2]);
                for ep_ix = 1:size(actv_epochs,1)
                    sig_chunk_ix = [find(hfa_actv.time==actv_epochs(ep_ix,1))...
                        find(hfa_actv.time==actv_epochs(ep_ix,2))];
                    % Report sign
                    if 0<=squeeze(mean(mean(hfa_actv.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
                        actv_ep_sign(ep_ix) = 1;
                    else
                        actv_ep_sign(ep_ix) = -1;
                    end
                end
                
                % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
                if any(actv_ep_sign==1)
                    %                 actv_code{roi_ix}    = [actv_code{roi_ix} 1];
                    actv_cnt(sbj_ix,roi_ix) = actv_cnt(sbj_ix,roi_ix)+1;
                end
                if any(actv_ep_sign==-1)
                    %                 actv_code{roi_ix}    = [actv_code{roi_ix} -1];
                    dact_cnt(sbj_ix,roi_ix) = dact_cnt(sbj_ix,roi_ix)+1;
                end
                both_ix = both_ix + 1;
            end
            
            % Check for RT correlations, get epochs
            if sum(squeeze(stat.mask(ch_ix,1,:)))>0
%                 mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
%                 mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
%                 % Convert to time
%                 for ep_ix = 1:size(mask_chunks,1)
%                     mask_chunks(ep_ix,1) = stat.time(mask_chunks(ep_ix,1));
%                     mask_chunks(ep_ix,2) = stat.time(mask_chunks(ep_ix,2));
%                 end
%                 %             cond_count(roi_ix) = cond_count(roi_ix) + 1;
%                 if strcmp('stim',event)
%                     % Exclude differences after the mean RT for this SBJ
%                     if min(mask_chunks(:,1)) < mean_RTs(sbj_ix)
%                         cond_g_count(sbj_ix,groi_ix) = cond_g_count(sbj_ix,groi_ix) + 1;
%                         both_ix = both_ix + 1;
%                     end
%                 else
                rt_cnt(sbj_ix,roi_ix) = rt_cnt(sbj_ix,roi_ix) + 1;
%                 end
            end
            
            % Check for ANOVA group effects
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(qvals(grp_ix,ch_ix,:))<0.05)
                    grp_cnt(sbj_ix,roi_ix,grp_ix) = grp_cnt(sbj_ix,roi_ix,grp_ix)+1;
                end
            end
        end
    end
    clear SBJ SBJ_vars w2 hfa_actv stat einfo actv_ch actv_ch_epochs tmp trl_info qvals
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
% Create and format the plot
fig_name = ['GRP_HFA_bar_perc_actv_' stat_id '_' roi_id '_' event_lab];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% % Create place holder line to initialize axes
% [ax,h1,h2] = plotyy(plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)],...
%        plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)]);
% delete(h1);
% delete(h2);

% Compile Data in Plotting Format
scat_vars = {actv_cnt, dact_cnt, rt_cnt};
bar_vars  = zeros([3+numel(grp_lab) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
    bar_vars(1,roi_ix) = sum(actv_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    bar_vars(2,roi_ix) = sum(dact_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    bar_vars(3,roi_ix) = sum(rt_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    for grp_ix = 1:numel(grp_lab)
        bar_vars(3+grp_ix,roi_ix) = sum(grp_cnt(:,roi_ix,grp_ix))/sum(elec_cnt(:,roi_ix));
        scat_vars{3+grp_ix} = squeeze(grp_cnt(:,:,grp_ix));
    end
end

% Plot Activations by ROI
b = {};
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(roi_list));   %bar for activation, deactivation, condition
for an_ix = 1:3+numel(grp_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{an_ix} = bar(bar_offsets+an_ix,diag(bar_vars(an_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(b{an_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        % Plot individual subject percentages as scatters on top
        has_elecs = squeeze(elec_cnt(:,roi_ix)~=0);
        if has_elecs
            s = scatter(scat_offsets(has_elecs)+bar_offsets(roi_ix)+an_ix,...
                scat_vars{an_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),50,'k*');
        end
    end
end
% if strcmp(event_lab,'stim')
%     legend([b{1},s],roi_list{:},'Individuals','Location','northwest');
% else
legend([b{1},s],roi_list{:},'Individuals','Location','northeast');
% end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = [0.5 3.5+numel(grp_lab)];
ax.XTick   = 1:3+numel(grp_lab);
ax.XColor  = 'k';
ax.XTickLabel = {'Activation','Deactivation','Corr(RT)',grp_lab{:}};

ax.YLabel.String   = '% Electrodes';
ax.YLabel.FontSize = 14;
% ax.YLim            = [0 ymaxs(plot_ix)];
ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Percentage of Electrodes Showing Significant Effects';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
ax.Title.FontSize = 16;

%% Save figure
if save_fig
    fig_dir = ['/home/knight/hoycw/PRJ_Error/results/HFA/GRP_summary_bar_ROI/actv_'...
        stat_id '_' roi_id '/' an_id '_mn' actv_win '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

end
