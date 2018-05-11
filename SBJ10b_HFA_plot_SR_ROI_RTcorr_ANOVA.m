function SBJ10b_HFA_plot_SR_ROI_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
                                                roi_id,plt_id,save_fig,fig_vis,fig_filetype)
% Plot time series for all significant effects by ROI
% clear all; %close all;
% fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Prep variables
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
if ~strcmp(plt_vars.sig_type,'bold')
    error('Only bold sig_type allowed');
end

% Get condition info
[grp_lab, ~, ~] = fn_group_label_styles(model_lab);
[rt_lab, ~, ~]     = fn_group_label_styles('RT');
cond_lab = {grp_lab{:}, rt_lab{:}};

% Get event timing
mean_RTs = zeros(size(SBJs));
event_lab = {'stim', 'resp'};

% Load all ROI info
[roi_list, roi_colors, einfo_roi_col] = fn_roi_label_styles(roi_id);

% Set up onset counts
ts       = {cell([numel(SBJs) numel(roi_list) numel(grp_lab)+1])...
            cell([numel(SBJs) numel(roi_list) numel(grp_lab)+1])};
sig_idx  = {cell([numel(SBJs) numel(roi_list) numel(grp_lab)+1])...
            cell([numel(SBJs) numel(roi_list) numel(grp_lab)+1])};
elec_roi_cnt  = zeros([2 numel(SBJs) numel(roi_list)]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load data
    %   !!! NEED TO CHANGE TO STAT_ID !!!
    f_name_s = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id_s '.mat'];
    f_name_r = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id_r '.mat'];
    tmp = load(f_name_s,'w2'); w2{1} = tmp.w2;
    tmp = load(f_name_r,'w2'); w2{2} = tmp.w2;
    tmp = load(f_name_s,'stat'); stat{1} = tmp.stat;
    tmp = load(f_name_r,'stat'); stat{2} = tmp.stat;
    clear tmp
    
    % Compute mean RT
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    mean_RTs(sbj_ix) = mean(trial_info.response_time);
    rt_plot_ix = find(stat{2}.time==0);
    
    sample_rate = (numel(stat{1}.time)-1)/(stat{1}.time(end)-stat{1}.time(1));
        
    %% Adjust Data
    % FDR correct ANOVA results
    qvals = {};
    win_lim = {}; win_center = {};
    for sr_ix = 1:numel(event_lab)
        % Convert % explained variance to 0-100 scale
        w2{sr_ix}.trial = w2{sr_ix}.trial*100;
        
        % FDR correct pvalues for ANOVA
        qvals{sr_ix} = NaN(size(w2{sr_ix}.pval));
        for ch_ix = 1:numel(stat{sr_ix}.label)
            [~, ~, ~, qvals{sr_ix}(:,ch_ix,:)] = fdr_bh(squeeze(w2{sr_ix}.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
        end
        
        % Get Sliding Window Parameters
        win_lim{sr_ix}    = fn_sliding_window_lim(stat{sr_ix}.time,win_len,win_step);
        win_center{sr_ix} = round(mean(win_lim{sr_ix},2));
    end
    
    % Trim data to plotting epoch
    %   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
    %   w2 should fit within that since it's averaging into a smaller window
    cfg_trim = [];
    cfg_trim.latency = plt_vars.plt_lim_S;
    stat{1} = ft_selectdata(cfg_trim,stat{1});
    cfg_trim.latency = plt_vars.plt_lim_R;
    stat{2} = ft_selectdata(cfg_trim,stat{2});
        
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
    if ~isempty(setdiff(stat{1}.label,einfo(:,1)))
        error('ERROR: Electrodes do not match between stat and einfo!');
    end
    
    %% Aggregate results per ROI
    for sr_ix = 1:numel(event_lab)
        for ch_ix = 1:numel(stat{sr_ix}.label)
            einfo_ix = strmatch(stat{sr_ix}.label(ch_ix),einfo(:,1),'exact');
            % Check for condition differences, get epochs
            if strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact')
                roi_ix = strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact');
                elec_roi_cnt(sr_ix,sbj_ix,roi_ix) = elec_roi_cnt(sr_ix,sbj_ix,roi_ix) + 1;
                
                % Get ANOVA group onsets
                for grp_ix = 1:numel(grp_lab)
                    if any(squeeze(qvals{sr_ix}(grp_ix,ch_ix,:))<0.05)
                        sig_idx{sr_ix}{sbj_ix,roi_ix,grp_ix} = [...
                                            sig_idx{sr_ix}{sbj_ix,roi_ix,grp_ix}; ...
                                            squeeze(qvals{sr_ix}(grp_ix,ch_ix,:)<0.05)'];
                        ts{sr_ix}{sbj_ix,roi_ix,grp_ix} = [...
                                            ts{sr_ix}{sbj_ix,roi_ix,grp_ix}; ...
                                            squeeze(w2{sr_ix}.trial(grp_ix,ch_ix,:))'];
                    end
                end
                
                % Get RT correlation onset
                if sum(squeeze(stat{sr_ix}.mask(ch_ix,1,:)))>0
                    sig_idx{sr_ix}{sbj_ix,roi_ix,numel(grp_lab)+1} = [...
                                            sig_idx{sr_ix}{sbj_ix,roi_ix,numel(grp_lab)+1}; ...
                                            squeeze(stat{sr_ix}.mask(ch_ix,:,:))'];
                    ts{sr_ix}{sbj_ix,roi_ix,numel(grp_lab)+1} = [...
                                            ts{sr_ix}{sbj_ix,roi_ix,numel(grp_lab)+1}; ...
                                            squeeze(stat{sr_ix}.rho(ch_ix,:,:))'];
                end
            end
        end
    end
    clear SBJ SBJ_vars hfa stat einfo w2
end

%% Plot GROI Results
% Find plot limits
ylims_grp = zeros([2 2]); yticks_grp = {};
for grp_ix = 1:numel(grp_lab)
    max_w2 = max([max(max(vertcat(ts{1}{:,:,grp_ix}))) ...
        max(max(vertcat(ts{2}{:,:,grp_ix})))]);
    min_w2 = min([min(min(vertcat(ts{1}{:,:,grp_ix}))) ...
        min(min(vertcat(ts{2}{:,:,grp_ix})))]);
    ylim_grp_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
    ylims_grp(grp_ix,:)  = [min_w2-ylim_grp_fudge max_w2+ylim_grp_fudge];
    yticks_grp{grp_ix} = 0:1:ylims_grp(grp_ix,2);
end

max_rho = max([max(max(vertcat(ts{1}{:,:,numel(grp_lab)+1}))) ...
               max(max(vertcat(ts{2}{:,:,numel(grp_lab)+1})))]);
min_rho = min([min(min(vertcat(ts{1}{:,:,numel(grp_lab)+1}))) ...
               min(min(vertcat(ts{2}{:,:,numel(grp_lab)+1})))]);
ylims_rt  = [round(min_rho*10)/10-0.1 round(max_rho*10)/10+0.1]; % extra on top and bottom for StdErr
yticks_rt = ylims_rt(1):0.1:ylims_rt(2);

% Create a figure for each channel
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    for cond_ix = 1:numel(cond_lab)
        grp_ix = strmatch(cond_lab{cond_ix},grp_lab);
        fig_name = [SBJ '_ANOVA_' stat_id '_SR_' roi_id '_' cond_lab{cond_ix}];
        figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 1 0.5],'Visible',fig_vis); %twice as wide for the double plot
        
        for sr_ix = 1:2
            subplot(1,2,sr_ix);
            hold on;
            roi_lines = cell(size(roi_list));
            roi_flag  = ones(size(roi_list));
            lgd_lab   = cell(size(roi_list));
            sig_elec_cnt = 0;
            
            % Plot effect size
            for roi_ix = 1:numel(roi_list)
                if strcmp(cond_lab{cond_ix},'RT') && ~isempty(ts{sr_ix}{sbj_ix,roi_ix,cond_ix})
                    tmp_line = plot(1:numel(ts{sr_ix}{sbj_ix,roi_ix,cond_ix}(1,:)),...
                        ts{sr_ix}{sbj_ix,roi_ix,cond_ix},...
                        'Color',roi_colors{roi_ix},'LineStyle','-');
                    
                    roi_lines{roi_ix} = tmp_line(1);
                    lgd_lab{roi_ix}   = [roi_list{roi_ix} '(n='...
                        num2str(size(ts{sr_ix}{sbj_ix,roi_ix,cond_ix},1)) '/'...
                        num2str(sum(elec_roi_cnt(sr_ix,sbj_ix,roi_ix))) ')'];
                    sig_elec_cnt = sig_elec_cnt + size(ts{sr_ix}{sbj_ix,roi_ix,cond_ix},1);
                elseif ~isempty(ts{sr_ix}{sbj_ix,roi_ix,cond_ix})
                    tmp_line = plot(win_center{sr_ix},...
                        ts{sr_ix}{sbj_ix,roi_ix,cond_ix},...
                        'Color',roi_colors{roi_ix},'LineStyle','-');
                    
                    roi_lines{roi_ix} = tmp_line(1);
                    lgd_lab{roi_ix}   = [roi_list{roi_ix} '(n='...
                        num2str(size(ts{sr_ix}{sbj_ix,roi_ix,cond_ix},1)) '/'...
                        num2str(sum(elec_roi_cnt(sr_ix,sbj_ix,roi_ix))) ')'];
                    sig_elec_cnt = sig_elec_cnt + size(ts{sr_ix}{sbj_ix,roi_ix,cond_ix},1);
                else
                    roi_flag(roi_ix) = 0;
                    lgd_lab{roi_ix}  = [];
                end
            end
%             lgd_lab(logical(roi_flag)) = [];
            
            % Plot Response Marker
            if strcmp(event_lab{sr_ix},'stim')
                x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
                event_time = -plt_vars.plt_lim_S(1)+mean_RTs(sbj_ix)*sample_rate;
                lgd_loc    = plt_vars.legend_loc_S;
            else
                x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
                event_time = rt_plot_ix;
                lgd_loc    = plt_vars.legend_loc_R;
            end
            if strcmp(cond_lab{cond_ix},'RT')
                event_line = line([event_time event_time],ylims_rt);
            else
                event_line = line([event_time event_time],ylims_grp(grp_ix,:));
            end
            set(event_line, 'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
            roi_lines = {roi_lines{:} event_line};
            roi_flag = [roi_flag 1];
            lgd_lab = {lgd_lab{:} 'RT'};
            
            % Plot significant time periods
            for roi_ix = 1:numel(roi_list)
                for ch_ix = 1:size(sig_idx{sr_ix}{sbj_ix,roi_ix,cond_ix},1)
                    if strcmp(cond_lab{cond_ix},'RT')
                        % RT significant epochs
                        sig_chunks = fn_find_chunks(squeeze(sig_idx{sr_ix}{sbj_ix,roi_ix,cond_ix}(ch_ix,:)));
                        sig_chunks(squeeze(sig_idx{sr_ix}{sbj_ix,roi_ix,cond_ix}(ch_ix,sig_chunks(:,1)))==0,:) = [];
                        for sig_ix = 1:size(sig_chunks,1)
                            sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
                            line(sig_times,squeeze(ts{sr_ix}{sbj_ix,roi_ix,cond_ix}(ch_ix,sig_times)),...
                                'Color',roi_colors{roi_ix},'LineStyle',plt_vars.sig_style,...
                                'LineWidth',plt_vars.sig_width);
                        end
                    else
                        % ANOVA factor significant periods
                        sig_chunks = fn_find_chunks(squeeze(sig_idx{sr_ix}{sbj_ix,roi_ix,cond_ix}(ch_ix,:)));
                        sig_chunks(squeeze(sig_idx{sr_ix}{sbj_ix,roi_ix,cond_ix}(ch_ix,sig_chunks(:,1)))==0,:) = [];
                        for sig_ix = 1:size(sig_chunks,1)
                            sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
                            line(win_center{sr_ix}(sig_times),...
                                squeeze(ts{sr_ix}{sbj_ix,roi_ix,cond_ix}(ch_ix,sig_times)),...
                                'Color',roi_colors{roi_ix},'LineStyle',plt_vars.sig_style,...
                                'LineWidth',plt_vars.sig_width);
                        end
                    end
                end
            end
            
            % Plotting parameters
            ax = gca;
            ax.Title.String  = [SBJ '-' cond_lab{cond_ix} ' (sig/total n='...
                num2str(sig_elec_cnt) '/' num2str(sum(elec_roi_cnt(sr_ix,sbj_ix,:))) '): ' event_lab{sr_ix}];
            if strcmp(cond_lab{cond_ix},'RT')
                ax.YLim          = ylims_rt;
                ax.YTick         = yticks_rt;
                ax.YLabel.String = 'Correlation with RT';
            else
                ax.YLim          = ylims_grp(grp_ix,:);
                ax.YTick         = yticks_grp{grp_ix};
                ax.YLabel.String = '% Variance Explained';
            end
            ax.XLim          = [0,size(sig_idx{sr_ix}{sbj_ix,1,numel(cond_lab)},2)];
            ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(sig_idx{sr_ix}{sbj_ix,1,numel(cond_lab)},2);
            ax.XTickLabel    = x_tick_lab;
            ax.XLabel.String = 'Time (s)';
            legend([roi_lines{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location',lgd_loc);
        end
        
        % Save figure
        if save_fig
            fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' stat_id '/SR/' an_id_s '-' an_id_r '/'];
            if ~exist(fig_dir,'dir')
                mkdir(fig_dir);
            end
            fig_filename = [fig_dir fig_name '.' fig_filetype];
            fprintf('Saving %s\n',fig_filename);
            saveas(gcf,fig_filename);
            %eval(['export_fig ' fig_filename]);
        end
    end
end

end
