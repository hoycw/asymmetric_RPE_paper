sbj_list = SBJs;
SBJ_colors = distinguishable_colors(numel(sbj_list));
cond_lab = {'easy','hard'};
acc_cond = zeros([numel(sbj_list) 2]);
for sbj_ix = 1:numel(sbj_list)
    SBJ = sbj_list{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'),'trl_info');
    for cond_ix = 1:2
        acc_cond(sbj_ix,cond_ix) = mean(trl_info.hit(strcmp(trl_info.cond,cond_lab{cond_ix})));
    end
end

%%
figure; hold on;
for sbj_ix = 1:numel(sbj_list)
    line([1 2],acc_cond(sbj_ix,:),'Color',SBJ_colors(sbj_ix,:));
    scatter([1 2], acc_cond(sbj_ix,:),30,SBJ_colors(sbj_ix,:));
end
ax = gca;
ax.XLabel.String = 'Condition';
ax.XTick = [1 2];
ax.XTickLabels = cond_lab;
ax.XLim = [0.5 2.5];

ax.YLim = [0 1];
ax.YLabel.String = 'Accuracy';

ax.Title.String = ['Patient Behavior (n = ' num2str(numel(sbj_list)) ')'];

set(ax,'FontSize',16');