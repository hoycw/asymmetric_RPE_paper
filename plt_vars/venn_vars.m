plt.fig_pos  = {[0 0 0.5 0.7], [0 0 0.5 0.8]};
plt.axis_vis = 'off';

% Venn features
plt.cmap       = parula(100);
plt.ang        = linspace(0, 2*pi);
plt.ax_fudge   = 3;
plt.alpha      = 0.6;
plt.line_color = [0 0 0];

% Text features
plt.title_sz     = 20;
plt.legend_sz    = 16;
plt.text_sz      = 20;
plt.text_align   = 'center';
plt.text_spacers = {...
        [-1 0; 1 0; ...%Main Effects
         0 -1],  % single pair
        [-1 0; 1 0; 0 1; ...%Main Effects
         0 -1; -1 1; 1 1; 0 0]};  % pairs + interaction
