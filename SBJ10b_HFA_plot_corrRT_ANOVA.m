function SBJ10b_HFA_plot_corrRT_ANOVA(SBJ,pipeline_id,stat_id,an_id,plt_id,save_fig,fig_vis,fig_filetype)
% Plots ANOVA results
% clear all; %close all;
% fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Error/scripts/');
addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
eval(['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Error/scripts/plt_vars/' plt_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');

f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
tmp = load(f_name,'w2'); w2 = tmp.w2;
% tmp = load(f_name_s,'hfa'); hfa = tmp.hfa;
tmp = load(f_name,'stat'); stat = tmp.stat;
clear tmp

% Load ROI and GM/WM info
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

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));

%% Prep Data
% FDR correct pvalues for ANOVA
win_lim = {}; win_center = {};
qvals = NaN(size(w2.pval));
for ch_ix = 1:numel(stat.label)
    pvals = squeeze(w2.pval(:,ch_ix,:));
    [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(pvals);%,0.05,'pdep','yes');
end

% Get Sliding Window Parameters
win_lim    = fn_sliding_window_lim(stat.time,win_len,win_step);
win_center = round(mean(win_lim,2));

% Convert % explained variance to 0-100 scale
w2.trial = w2.trial*100;

% Trim data to plotting epoch
%   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
%   w2 should fit within that since it's averaging into a smaller window
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim;
% hfa{1}  = ft_selectdata(cfg_trim,hfa{1});
stat = ft_selectdata(cfg_trim,stat);

% % Compute mean RT per condition
% RTs = round(1000*trl_info.rt); % converts sec to ms
% for cond_ix = 1:numel(grp_lab{congr_ix})
%     RT_means{cond_ix} = mean(RTs(fn_condition_index([grp_lab{congr_ix}{cond_ix}], trl_info.condition_n)==1));
%     % Add in the baseline offset to plot correctly
%     RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
% end

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Error/results/HFA/' SBJ '/' stat_id '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Find plot limits
max_w2 = max(max(max(w2.trial)));
min_w2 = min(min(min(w2.trial)));
ylim1_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims1  = [min_w2-ylim1_fudge max_w2+ylim1_fudge];
yticks1 = 0:1:ylims1(2);

max_rho = max(max(squeeze(stat.rho)));
min_rho = min(min(squeeze(stat.rho)));
ylims2  = [round(min_rho*10)/10-0.1 round(max_rho*10)/10+0.1]; % extra on top and bottom for StdErr
yticks2 = ylims2(1):0.1:ylims2(2);
y_sig = zeros([1 numel(grp_lab)+1]);
y_sig(1) = mean([min_w2,max_w2]);
for grp_ix = 2:numel(grp_lab)+2
    y_sig(grp_ix) = y_sig(grp_ix-1)+ylim1_fudge;
end

% Create a figure for each channel
sig_ch = zeros([numel(stat.label),1+numel(grp_lab)]);
for ch_ix = 1:numel(stat.label)
    fig_name = [SBJ '_' stat_id '_' stat.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.8],'Visible',fig_vis); %twice as wide for the double plot
    
    hold on;
    main_lines = [];
    
    % Plot var_exp
    [axs,grp_lines,rt_line] = plotyy(...
        repmat(win_center,size(grp_lab)),squeeze(w2.trial(:,ch_ix,:))',...
        1:numel(stat.time),squeeze(stat.rho(ch_ix,:,:))');
    for grp_ix = 1:numel(grp_lab)
        set(grp_lines(grp_ix),'Color',[grp_colors{grp_ix}],'LineStyle',grp_style{grp_ix});
    end
    set(rt_line,'Color',rt_color{1},'LineStyle',rt_style{1});
    main_lines = [main_lines grp_lines' rt_line];
    lgd_lab = {grp_lab{:} 'corr(RT)'};
    
    % Plot events: stim, target, feedback onset, feedback offset
    x_tick_lab       = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
%     event_lab    = {'stim','target','fb on','fb off'};
    event_styles = {'--', '-', '-'};
    event_times = [(trl_info.timing.target-plt_vars.plt_lim(1))*sample_rate...
        (trl_info.timing.target+trl_info.timing.fb_delay-plt_vars.plt_lim(1))*sample_rate...
        (trl_info.timing.trl_len-plt_vars.plt_lim(1))*sample_rate];
    event_lines = [];
    for event_ix = 1:numel(event_times)
        event_lines(event_ix) = line([event_times(event_ix) event_times(event_ix)],ylims1,...
            'LineWidth',plt_vars.evnt_width,'Color','k','LineStyle',event_styles{event_ix});
    end
%         % Plot Response Marker
%         event_line = line([find(stat.time==0) find(stat.time==0)],ylims1,...
%             'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%             'LineStyle',plt_vars.evnt_style);
%         main_lines = [main_lines event_line];
    
    % Plot significant time periods
    for grp_ix = 1:numel(grp_lab)+1
        if grp_ix <= numel(grp_lab)
            if any(squeeze(qvals(grp_ix,ch_ix,:))<0.05)
                % Track channels with significant var_exp
                sig_ch(ch_ix,grp_ix) = 1;
                
                % Find significant periods
                if strcmp(plt_vars.sig_type,'bold')
                    sig_chunks = fn_find_chunks(squeeze(qvals(grp_ix,ch_ix,:))<0.05);
                    sig_chunks(squeeze(qvals(grp_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
                    set(gcf,'CurrentAxes',axs(1));
                    for sig_ix = 1:size(sig_chunks,1)
                        if diff(sig_chunks(sig_ix,:))==0
                            scatter(win_center(sig_chunks(sig_ix,1)),squeeze(w2.trial(grp_ix,ch_ix,sig_chunks(sig_ix,1))),...
                                plt_vars.sig_scat_size,grp_colors{grp_ix},plt_vars.sig_scat_mrkr,'filled');
                        else
                            line(win_center(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                                squeeze(w2.trial(grp_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                                'Color',[grp_colors{grp_ix}],'LineStyle',plt_vars.sig_style,...
                                'LineWidth',plt_vars.sig_width);
                        end
                    end
                elseif strcmp(plt_vars.sig_type,'scatter')
                    sig_times = win_center(squeeze(qvals(grp_ix,ch_ix,:))<0.05);
                    scatter(sig_times,repmat(y_sig(grp_ix),size(sig_times)),...
                        plt_vars.sig_scat_size,grp_colors{grp_ix},plt_vars.sig_scat_mrkr);
                elseif strcmp(plt_vars.sig_type,'patch')
                    sig_chunks = fn_find_chunks(squeeze(qvals(grp_ix,ch_ix,:))<0.05);
                    sig_chunks(squeeze(qvals(grp_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
                    error('sig_type = patch needs sig_times variable!');
                    % Plot Significance Shading
                    fprintf('%s %s -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                        stat.label{ch_ix},grp_lab{grp_ix},size(sig_chunks,1));
                    for sig_ix = 1:size(sig_chunks,1)
                        sig_times = win_lim(sig_chunks(sig_ix,:),:);
                        patch([sig_times(1,1) sig_times(1,1) sig_times(2,2) sig_times(2,2)], ...
                            [ylims1(1) ylims1(2) ylims1(2) ylims1(1)],...
                            grp_colors{grp_ix},'FaceAlpha',plt_vars.sig_alpha);
                    end
                end
            else
                fprintf('%s %s -- NO SIGNIFICANT CLUSTERS FOUND...\n',...
                    stat.label{ch_ix},grp_lab{grp_ix});
            end
        else
            % RT correlation significant time periods
            if sum(stat.mask(ch_ix,:))>0
                % Track channels with significant var_exp
                sig_ch(ch_ix,grp_ix) = 1;
                
                sig_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,:,:)));
                sig_chunks(squeeze(stat.mask(ch_ix,:,sig_chunks(:,1)))==0,:) = [];
                sig_times = find(squeeze(stat.mask(ch_ix,:,:)==1));
                set(gcf,'CurrentAxes',axs(2));
                if strcmp(plt_vars.sig_type,'bold')
                    for sig_ix = 1:size(sig_chunks,1)
                        sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
                        line(sig_times,squeeze(stat.rho(ch_ix,:,sig_times)),...
                            'Color',rt_color{1},'LineStyle',plt_vars.sig_style,...
                            'LineWidth',plt_vars.sig_width);
                    end
                elseif strcmp(plt_vars.sig_type,'scatter')
                    scatter(sig_times,repmat(y_sig(grp_ix),size(sig_times)),...
                        plt_vars.sig_scat_size2,rt_color{1},plt_vars.sig_scat_mrkr2);
                end
                fprintf('%s RT -- SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                    stat.label{ch_ix});%,size(sig_chunks,1));
            else
                fprintf('%s RT -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',...
                    stat.label{ch_ix});
            end
        end
    end
    
    % Plotting parameters
    axs(1).Title.String  = strcat(stat.label{ch_ix}, ' (', einfo(ch_ix,2), ')');
    axs(1).Box           = 'off';
    axs(1).YLim          = ylims1;
    axs(1).YTick         = yticks1;
    axs(1).YLabel.String = '% Variance Explained';
    axs(1).XLim          = [0,size(stat.time,2)];
    axs(1).XTick         = 0:plt_vars.x_step_sz*sample_rate:size(stat.time,2);
    axs(1).XTickLabel    = x_tick_lab;
    axs(1).XLabel.String = 'Time (s)';
    axs(2).YColor        = 'k';
    axs(2).YLim          = ylims2;
    axs(2).YTick         = yticks2;
    axs(2).YLabel.String = 'Correlation with RT';
    legend(main_lines,lgd_lab{:},'Location',plt_vars.legend_loc);
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

%% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'ch_sig_list.csv'];
sig_report = fopen(sig_report_filename,'w');
fprintf(sig_report,'%s: %s for %s \n',SBJ,stat_id,an_id);
fprintf(sig_report,'label,%s,RT\n',strjoin(grp_lab,','));
for ch_ix = 1:numel(stat.label)
    fprintf(sig_report,'%s,',stat.label{ch_ix});
    for grp_ix = 1:numel(grp_lab)
        fprintf(sig_report,'%.0f,',sig_ch(ch_ix,grp_ix));
    end
    fprintf(sig_report,'%.0f\n',sig_ch(ch_ix,grp_ix+1));
end
fclose(sig_report);

end
