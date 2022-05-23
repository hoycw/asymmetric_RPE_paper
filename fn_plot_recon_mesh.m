function fig = fn_plot_recon_mesh(elec, roi_mesh, roi_mesh_lab, rcn, plot_name)
%% Plot 3D surface mesh reconstructions with colored electrodes
% INPUTS:
%   elec [struct] - Fieldtrip elec struct with locations and colors to plot
%   roi_mesh [cell array] - struct mesh(es) to plot (.pos and .tri fields)
%   roi_mesh_lab [cell array] - string labels for each mesh
%   rcn [struct] - recon varaibles to control plotting
%   plot_name [str] - name of the figure
% OUTPUTS:
%   fig [struct] - figure handle for the plot

fig = figure('Name',plot_name);

% Plot 3D mesh(es)
for mesh_ix = 1:length(roi_mesh)
    fprintf('Plotting %s mesh (%i/%i) for %s...\n',roi_mesh_lab{mesh_ix},mesh_ix,length(roi_mesh),rcn.plot_roi);
    ft_plot_mesh(roi_mesh{mesh_ix}, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', rcn.mesh_alpha);
end

% Plot electrodes on top (one by one to get different colors)
for e = 1:numel(elec.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.color, 'label', rcn.lab_arg);
end

% Set viewing properties
view(rcn.view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(fig, 'windowkeypressfcn',   @cb_keyboard);

end