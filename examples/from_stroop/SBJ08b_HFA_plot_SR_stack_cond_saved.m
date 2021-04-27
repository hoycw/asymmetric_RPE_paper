function SBJ08b_HFA_plot_SR_stack_cond_saved(SBJ,conditions,an_id_s,an_id_r,actv_win,...
                                        plt_id,save_fig,fig_vis,fig_filetype)
% Plots single trial stack for both stimulus- and response-locked HFA computed in SBJ08a_HFA_actv
%   sorts by condition, then by RT; scatter for RTs in stim-locked
% clear all; %close all;
!!! not converted from PRJ_Stroop!
if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

% Load data
hfa_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id_s,'.mat');
hfa_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id_r,'.mat');
tmp = load(hfa_filename1,'hfa'); hfa{1} = tmp.hfa;
tmp = load(hfa_filename2,'hfa'); hfa{2} = tmp.hfa;
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));
if ~isempty(setdiff(hfa{1}.label,hfa{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2} = ft_selectdata(cfg_trim,hfa{2});

% Compile cond_type, RT, trial_n
[cond_lab, cond_colors, ~] = fn_condition_label_styles(conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end
event_lab = {'stim', 'resp'};

[cond_mat,~] = find(cond_idx); 
cond_mat = horzcat(cond_mat,round(1000*trial_info.response_time),[1:numel(cond_mat)]');
cond_mat = sortrows(cond_mat,[1 2]);
RT_mean = mean(round(1000*trial_info.response_time)); % converts sec to ms
% Add in the baseline offset to plot correctly
RT_mean = RT_mean-plt_vars.plt_lim_S(1)*1000;

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/stack_' conditions '/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa{1}.label)
    % Plot parameters
    fig_name = [SBJ '_' conditions '_SR_stack_' hfa{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
    
    % Get color limits
    clims = NaN([2 2]);
    for sr_ix = 1:2
        clims(sr_ix,1) = prctile(reshape(hfa{sr_ix}.powspctrm(:,ch_ix,:,:),...
            [1 size(cond_mat,1)*numel(hfa{sr_ix}.time)]),plt_vars.clim_perc(1));
        clims(sr_ix,2) = prctile(reshape(hfa{sr_ix}.powspctrm(:,ch_ix,:,:),...
            [1 size(cond_mat,1)*numel(hfa{sr_ix}.time)]),plt_vars.clim_perc(2));
    end
    clims = [min(clims(:,1)) max(clims(:,2))];
        
    for sr_ix = 1:numel(event_lab)
        subplot(1,numel(event_lab),sr_ix);
        hold on;
        ax     = gca;
        
        % Plot Single Trials Per Condition
        imagesc(squeeze(hfa{sr_ix}.powspctrm(cond_mat(:,3),ch_ix,:,:)));
        set(gca,'YDir','normal');
        if strcmp(event_lab{sr_ix},'stim')
            x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            event_time = -plt_vars.plt_lim_S(1)*sample_rate;
            for cond_ix = 1:numel(cond_lab)
                idx = cond_mat(:,1)==cond_ix;
                scat = scatter(cond_mat(idx,2)-plt_vars.plt_lim_S(1)*sample_rate,find(idx),...
                    'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k');
            end
        else
            x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            event_time = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        ylim([1 size(cond_mat,1)]);
        event_line = line([event_time event_time],ylim,...
            'LineWidth',plt_vars.evnt_width,'Color','k');
        
        % Plotting parameters
        ax = gca;
%         ax.legend = plt_vars.legend;
        ax.Title.String  = [hfa{sr_ix}.label{ch_ix} ': ' event_lab{sr_ix} ' trials'];
        ax.XLim          = [0,numel(hfa{sr_ix}.time)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:numel(hfa{sr_ix}.time);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Trials';
%         legend([roi_lines{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location',lgd_loc);
        cbar = colorbar;
        caxis(clims);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
