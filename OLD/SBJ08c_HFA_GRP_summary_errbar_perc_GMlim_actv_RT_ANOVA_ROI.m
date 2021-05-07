function SBJ08c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,stat_id,proc_id,an_id,actv_win,roi_id,...
                                                            atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for active, RT correlation, and ANOVA epochs
%   RT correlation: any significance in stat_lim
%   ANOVA factors: any significance in stat_lim, after FDR correction
% OUTPUTS:
%   Bar chart with % active, % deactivated, % RT correlations, % ANOVA factors
% clear all; %close all;
% fig_ftype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
if ischar(plot_scat); plot_scat = str2num(plot_scat); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Prep variables
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
% plt_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/plt_vars/' plt_id '_vars.m'];
% eval(plt_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
% if rt_correlation
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
% end

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);

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
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');
    if strcmp(event_type,'stim')
        mean_RTs(sbj_ix) = mean(trl_info.rt);
    end
    % Load ANOVA
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',stat_id,'_',an_id,'.mat'));
    
    % Load actv stats
    actv_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_actv_mn',actv_win,'.mat');
    load(actv_fname,'actv_ch','actv_ch_epochs');
    
    % Load HFA (to get sign of activation)
    hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
    tmp = load(hfa_fname,'hfa'); hfa = tmp.hfa;
    
    %% Load ROI and GM/WM info
    if strcmp(atlas_id,'Yeo7')
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_mni_v_' atlas_id '.mat'];
    else
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat'];
    end
    load(elec_fname);
    
    % HACK!!!
    if strcmp(SBJ,'CP24')
        cfgs = [];
        cfgs.channel = {'all','-RTO4'};
        stat = ft_selectdata(cfgs,stat);
        hfa = ft_selectdata(cfgs,hfa);
    elseif strcmp(SBJ,'IR57')
        cfgs = [];
        cfgs.channel = {'all','-LAM4-5','-LAM5-6','-RIN8-9','-RIN9-10',...
                        '-RSM1-2','-RSM2-3','-RSM3-4','-RTI9-10'};
        stat = ft_selectdata(cfgs,stat);
        hfa = ft_selectdata(cfgs,hfa);
    elseif strcmp(SBJ,'IR68')
        cfgs = [];
        cfgs.channel = {'all','-LPC5-6','-LPC6-7','-LPC7-8'};
        stat = ft_selectdata(cfgs,stat);
        hfa = ft_selectdata(cfgs,hfa);
    end

    % Sort elecs by stat labels
    if ~strcmp(SBJ,'IR66')  % HACK!!!
        cfgs = []; cfgs.channel = stat.label;
        elec = fn_select_elec(cfgs,elec);
        elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    else
        elec.roi = elec.(roi_id);
    end
    
    % Exclude elecs not in atlas ROIs
    if ~plot_out
        cfgs = []; cfgs.channel = fn_select_elec_lab_match(elec, 'b', atlas_id, roi_id);
        elec = fn_select_elec(cfgs, elec);
        stat = ft_selectdata(cfgs,stat);
        w2   = ft_selectdata(cfgs,w2);
    end
    
    % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
    gm_bin = ones(size(elec.label));
%     if gm_thresh>0
%         gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
%     else
%     end
    
    %% Process parameters
    %!!! is this the best way to do this??? Maybe not...
    sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
    
    % Restrict hfa to stat_lim (should obviate some of below comparisons)
    cfg_lim = [];
    cfg_lim.latency = stat_lim;
    hfa      = ft_selectdata(cfg_lim,hfa);
    stat     = ft_selectdata(cfg_lim,stat);
    
    % FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    % Confirm channels and time axis are the same
    %   All the rounding and such is because some stupid rounding errors...
    same_start = round(hfa.time(1)*uint8(sample_rate))==round(stat.time(1)*uint8(sample_rate));
    same_end   = round(hfa.time(end)*uint8(sample_rate))==round(stat.time(end)*uint8(sample_rate));
    same_numel = size(hfa.time,2)==size(stat.time,2);
    if ~same_start || ~same_end || ~same_numel
        error('time axes are not the same across hfa analyses!');
    end
    if ~isempty(setdiff(stat.label,hfa.label));
        error('Different electrodes across hfa analyses!');
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        % If elec matches roi_list, get stats
        if any(strcmp(elec.roi{ch_ix},roi_list)) && gm_bin(ch_ix)
            roi_ix = find(strcmp(elec.roi{ch_ix},roi_list));
            elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
            
            % Check if active, get epochs
            if any(strcmp(stat.label{ch_ix},actv_ch))
                % Find significant epoch indices
                actv_epochs = actv_ch_epochs{strcmp(stat.label{ch_ix},actv_ch)};
                
%                 % Toss late epochs if S-locked
%                 if strcmp('stim',event_type)
%                     actv_epochs(actv_epochs(:,1)>mean_RTs(sbj_ix),:) = [];
%                 end
                
                % Find sign of (de)activation
                actv_ep_sign = NaN([1 size(actv_epochs,1)]);
                sig_chunk_ix = NaN([1 2]);
                for ep_ix = 1:size(actv_epochs,1)
                    sig_chunk_ix = [find(hfa.time==actv_epochs(ep_ix,1))...
                        find(hfa.time==actv_epochs(ep_ix,2))];
                    % Report sign
                    if 0<=squeeze(mean(mean(hfa.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
                        actv_ep_sign(ep_ix) = 1;
                    else
                        actv_ep_sign(ep_ix) = -1;
                    end
                end
                
                % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
                if any(actv_ep_sign==1)
                    actv_cnt(sbj_ix,roi_ix) = actv_cnt(sbj_ix,roi_ix)+1;
                end
                if any(actv_ep_sign==-1)
                    dact_cnt(sbj_ix,roi_ix) = dact_cnt(sbj_ix,roi_ix)+1;
                end
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
%                 if strcmp('stim',event_type)
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
    clear SBJ SBJ_vars w2 stat elec hfa actv_ch actv_ch_epochs tmp trl_info qvals
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
if plot_scat
    scat_suffix = '_SBJscat';
else
    scat_suffix = '';
end
% Create and format the plot
fig_name = ['GRP_HFA_errbar_perc_GMlim' num2str(gm_thresh*100) '_actv_' stat_id '_' roi_id '_' event_type scat_suffix];
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
bar_data  = zeros([3+numel(grp_lab) numel(roi_list)]);
var_data  = zeros([3+numel(grp_lab) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
%     bar_vars(1,roi_ix) = sum(cI_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    bar_data(1,roi_ix) = sum(actv_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(1,roi_ix) = nanstd(actv_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    bar_data(2,roi_ix) = sum(dact_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(2,roi_ix) = nanstd(dact_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    bar_data(3,roi_ix) = sum(rt_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(3,roi_ix) = nanstd(rt_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    for grp_ix = 1:numel(grp_lab)
        bar_data(3+grp_ix,roi_ix) = sum(grp_cnt(:,roi_ix,grp_ix))/sum(elec_cnt(:,roi_ix));
        var_data(3+grp_ix,roi_ix) = nanstd(grp_cnt(:,roi_ix,grp_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
        scat_vars{3+grp_ix} = squeeze(grp_cnt(:,:,grp_ix));
    end
end

% Plot Activations by ROI
b = {};
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(roi_list));   %bar for activation, deactivation, condition
for cond_ix = 1:3+numel(grp_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{cond_ix} = bar(bar_offsets+cond_ix,diag(bar_data(cond_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(b{cond_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        % Plot individual subject percentages as scatters on top
        line([bar_offsets(roi_ix) bar_offsets(roi_ix)]+cond_ix,...
            [bar_data(cond_ix,roi_ix)+var_data(cond_ix,roi_ix) bar_data(cond_ix,roi_ix)-var_data(cond_ix,roi_ix)],...
            'Color','k','LineWidth',1.5);
        % Overlay scatter for individual SBJ
        if plot_scat
            has_elecs = squeeze(elec_cnt(:,roi_ix)~=0);
            s = scatter(scat_offsets(has_elecs)+bar_offsets(roi_ix)+cond_ix,...
                scat_vars{cond_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),50,'k*');
        end
    end
end
% if strcmp(event_type,'stim')
%     leg_loc = 'northwest';
% else
%     leg_loc = 'northeast';
% end
leg_loc = 'northeast';
if plot_scat
    legend([b{1},s],roi_list{:},'Individuals','Location',leg_loc);
else
    legend(b{1},roi_list{:},'Location',leg_loc);
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = [0.5 3.5+numel(grp_lab)];
ax.XTick   = 1:3+numel(grp_lab);
ax.XColor  = 'k';
ax.XTickLabel = [{'Activation','Deactivation','Corr(RT)'},grp_lab];

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLabel.FontSize = 14;
ax.YLim            = [0 ax.YLim(2)];%[0 0.6];%
ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Proportion of Electrodes Showing Significant Effects';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
ax.Title.FontSize = 16;

%% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP_summary_errbar_ROI/actv_'...
        stat_id '_' roi_id '/' an_id '_mn' actv_win '/'];
    if ~exist(fig_dir,'dir')
        [~,~] = mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
