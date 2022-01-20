if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);

%% MI lag 0 on HFA
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';
roi_id = 'MPFCINS';
atlas_id = 'Dx';
chansel = 'all';%active';%'active';%'maxactive';
method = 'MI';
Nbins = 2.^[2:1:6];
for s = 1%2:numel(SBJs)
  SBJ = SBJs{s};
  SBJ10a_conn_lag0(SBJ, method, proc_id, an_id, roi_id, atlas_id, Nbins, chansel)
end
%% corr lag 0 on HFA
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';
roi_id = 'MPFCINS';
atlas_id = 'Dx';
chansel = 'all';%active';%'active';%'maxactive';
method = 'corr';
Nbins = 1;%2.^[2:1:6];
for s = 1%1:numel(SBJs)
  SBJ = SBJs{s};
  SBJ10a_conn_lag0(SBJ, method, proc_id, an_id, roi_id, atlas_id, Nbins, chansel)
end
%% Xcorr on HFA
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';
roi_id = 'MPFCINS';
atlas_id = 'Dx';
maxlag = 0.6; % in seconds
chansel = 'all';%'active';%'maxactive';
stepsize = 0.05;
for s = 1%:numel(SBJs)
  SBJ = SBJs{s};
  SBJ10b_conn_Xcorr(SBJ, proc_id, an_id, roi_id, atlas_id, maxlag, stepsize, chansel)
end
%% Windowed MI on HFA
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';
roi_id = 'MPFCINS';
atlas_id = 'Dx';
stepsize = 0.05; % in seconds
tlags = [0,0.6]; % lag limits in seconds;
twin = [0,0.6]; % time window
Nbins = 8; %2.^[2:1:6];
chansel = 'all';%'active';%'maxactive';
method = 'MI';%'MI';
for s = 1%4:numel(SBJs)
  SBJ = SBJs{s};
  SBJ10c_conn_win(SBJ, method, proc_id, an_id, roi_id, atlas_id, twin, tlags, stepsize, Nbins, chansel)
end
%% Windowed corr on HFA
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';
roi_id = 'MPFCINS';
atlas_id = 'Dx';
stepsize = 0.05; % in seconds
tlags = [0,0.6]; % lag limits in seconds;
twin = [0,0.6]; % time window
Nbins = []; %2.^[2:1:6];
chansel = 'all';%'active';%'maxactive';
method = 'corr';%'MI';
for s = 1%2:numel(SBJs)
  SBJ = SBJs{s};
  SBJ10c_conn_win(SBJ, method, proc_id, an_id, roi_id, atlas_id, twin, tlags, stepsize, Nbins, chansel)
end