plt.plt_lim    = [0 1];
% plt.x_data     = [0:0.025:1.5];
% plt.x_step_sz  = 0.2;
plt.legend     = 0;
plt.legend_loc = 'northeast';

%plt.actv_color  = {'r','b'}; % negative, positive
%plt.actv_alpha  = 0.5;
%plt.actv_height = 0.5;

plt.cond_color  = 'k';
plt.cond_width  = 2.5;

plt.grp_metric         = 'all';
% Only set this to override the default (ROI or reg_colors)
%plt.violin_scat_colors = 'gROI';

plt.roi_width = 2;
plt.roi_color = 'k';
plt.roi_style = ':';

plt.groi_width = 4;
plt.groi_color = 'k';
plt.groi_style = ':';
