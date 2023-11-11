function SBJ11f_HFA_conn_plot_grp_mLME_recon(SBJ_id, proc_id, an_id, model_id, conn_id,stat_id,...
    atlas_id, roi_id, rcn, swap_Xcorr, varargin)
%% Plot a reconstruction with electrodes colored according to statistics
%   This version uses categories based on EpnRPE model
% INPUTS:
%   SBJ_id [str] - ID of subject list to load
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   model_id [str] - ID of the model to load
%   stat_id [str] - ID of the stats (e.g., 'mLME_St0t6_WL05_WS25')
%   cat_id [str] - ID of the categories to group the coefficients
%       should be 'puns' for pos, neg, signed, unsigned/salience
%   atlas_id [str] - ID of atlas to select ROIs: {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%   rcn [struct] - options for recon plotting, see fn_process_recon_vars.m
%       .plot_roi [str] - which surface mesh to plot ('INS','MPFC',etc.)
%       .hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%       .mirror [0/1] - mirror elecs from one hemisphere to the other
%   OPTIONAL:
%       skip_reg [str] - name of one regressor to zero out and skip (not plot)

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Process plotting params
if ~isstruct(rcn); error('rcn is not a struct!'); end

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            rcn.view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            rcn.mesh_alpha = varargin{v+1};
        elseif strcmp(varargin{v},'skip_reg') && ischar(varargin{v+1})
            skip_reg = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype')
            fig_ftype = varargin{v+1};
        elseif strcmp(varargin{v},'fig_res')
            fig_res = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Implement the default options
if ~exist('save_fig','var');    save_fig = 1; end
if ~exist('fig_ftype','var');   fig_ftype = 'fig'; end
if ~exist('fig_res','var');     fig_res = 300; end
% ROI info
rcn = fn_process_recon_vars(rcn);
% rcn.view_angle = [-90,0];
%rcn.view_angle = [-130,45];

[~, ~, roi_field] = fn_roi_label_styles(roi_id);

%% Load data
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
%eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
if strcmp(mdl.model_lab,'RL3D')
    [reg_lab, reg_names, reg_colors, reg_styles] = fn_performanceRL_regressor_label_styles(mdl.model_lab);
else
    [reg_lab, reg_names, reg_colors, reg_styles, reg_mrkrs] = fn_regressor_label_styles(mdl.model_lab);
end
[~, names, colors, ~, ~] = fn_puns_category_label_styles('puns');

reg_cmaps = {'gray','autum','winter'};


if numel(reg_lab) < 2 || numel(reg_lab) > 3; error('why venn?'); end
%[cat_lab, cat_names, cat_colors, ~, ~] = fn_puns_category_label_styles(cat_id);

SBJs = fn_load_SBJ_list(SBJ_id);

load([root_dir 'PRJ_Error/data/GRP/stats/main_ft_' model_id '_' an_id '_' stat_id '_' conn_id  '_chancoef.mat']);

% Identify analysis type
if contains(an_id,'ERP')
    an_dir = 'ERP';
elseif contains(an_id,'HG')
    an_dir = 'HFA';
elseif  contains(an_id,'TFR')
    an_dir = 'TFR';
else
    error('unknown an_id');
end

if strcmp(rcn.plot_roi,'INS')
    roi_label = 'INS';
elseif strcmp(rcn.plot_roi,'MPFC')
    roi_label = 'MPFC';
elseif strcmp(rcn.plot_roi,'MPFCINS')
    roi_label = {'MPFC','INS'};
else
    error('SBJ08g not run yet for anything except MPFCINS');
end
%%
%% Get significant channel labels per regressor, peaks and latencies
regs = conn_stats_chan.feature;
% sig_reg_elecs = cell(size(regs));
% peaks = cell(size(regs));
% lats = cell(size(regs));
% dirs = cell(size(regs));
% for reg_ix = 1:length(regs)
%     elec_ix = find(sum(conn_stats_chan.pvals{1}(:,reg_ix+1,:) < .05,3) > 0);
%     % find peaks, lats and labels
%     peaks{reg_ix} = NaN(size(elec_ix));
%     lats{reg_ix} = NaN(size(elec_ix));
%     dirs{reg_ix} = NaN(size(elec_ix));
%     sig_reg_elecs{reg_ix} = cell(size(elec_ix,1),2);
%     for ch =  1:length(elec_ix)
%         %select channel and regressor data:
%         cchdata = squeeze(conn_stats_chan.coefs{1}(elec_ix(ch),reg_ix+1,:));
%         
%         %finde peaks latencies and directionality
%         [cpeak, latix] = findpeaks(abs(cchdata));
%         [~,pidx] = max(cpeak);
%         peaks{reg_ix}(ch) = cchdata(latix(pidx(1)));
%         lats{reg_ix}(ch) =  conn_stats_chan.time(latix(pidx(1)));
%         dirs{reg_ix}(ch) = sign(lats{reg_ix}(ch));
%         
%         % obtain channel labels
%         cur_lab = conn_stats_chan.chan_label{1}{elec_ix(ch)};
%         lab_split = strsplit(cur_lab,{' ','_to_'});
%         sig_reg_elecs{reg_ix}{ch,1} = [lab_split{1} ' ' lab_split{2}];
%         sig_reg_elecs{reg_ix}{ch,2} = [lab_split{1} ' ' lab_split{3}];
%     end
% end
%% Load elec and select within ROI and hemisphere
[elec_sbj, ~] = fn_load_grp_elec_ROI(SBJs,proc_id,atlas_id,roi_id,rcn);

%% Combine elec structs
elec = ft_appendsens([],elec_sbj{:});

%% Create 3D mesh based on atlas ROIs
[roi_mesh, roi_mesh_lab] = fn_load_recon_mesh_ROI(atlas_id,rcn.plot_roi,rcn.hemi);

%% align electrodes with mesh
elec = fn_align_elec2mesh(elec,roi_mesh);
%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
save_fig = true;
if save_fig
    out_dir = [root_dir 'PRJ_Error/results/' an_dir '/GRP/recon_conn/' model_id '/' an_id '/' stat_id '/' conn_id '/'];
    if ~exist(out_dir,'dir'); [~] = mkdir(out_dir); end
end
cat_labels = {'pRPE','nRPE','sRPE','uRPE'};
%cat_maps = {'autumn','winter','copper','summer'};
cat_maps = {{[255,0,0],[255,128,0]},...
            {[0,0,204],[0,204,204]},...
            {[209,151,105],[178,32,102]},...
            {[0,158,115],[153,153,0]}};
map_dir = [-1,-1,1,-1];
regs = {'p','n'};
reg_names = {'positive RPE','negative RPE'};
ctable = conn_stats_chan.peaks{1,1};

% split electrode names
split_elec = cellfun(@(x) strsplit(x,{' ','_to_'}),ctable.pair,'UniformOutput',false);
ctable.elec_from = cellfun(@(x) [x{1},' ',x{2}],split_elec,'UniformOutput',false);
ctable.elec_to = cellfun(@(x) [x{1},' ',x{3}],split_elec,'UniformOutput',false);

for cat = 1:length(cat_labels)
    ccat = cat_labels{cat};
    cmap = cat_maps{cat};
    cmap_dir = map_dir(cat);
    cat_name = names{cat};
    for reg = 1:length(regs)
        creg = regs{reg};
        creg_name = reg_names{reg};
        dirs = sign(ctable.([creg 'lat']));
        udirs = unique(dirs);
        disp(udirs)
        for dd = 1:length(udirs)
            celec = elec;
            d = udirs(dd);
            if swap_Xcorr == 1
                d = d*-1;
            end
            cdir_ix = ismember(dirs,d) & ismember(ctable.category,ccat);
            cur_reg_elecs = [ctable.elec_from(cdir_ix),ctable.elec_to(cdir_ix)];
            sig_ix = find(ismember(celec.label, cur_reg_elecs));
            celec.chanpos = celec.chanpos(sig_ix,:);
            celec.chantype = celec.chantype(sig_ix,:);
            celec.chanunit = celec.chanunit(sig_ix,:);
            celec.label= celec.label(sig_ix,:);
            celec.color = repmat(colors{cat},length(sig_ix),1);
            %celec.color = repmat([217 217 214]/255, size(celec.chanpos,1), 1);
            %celec.color(sig_ix,:) = repmat(reg_colors{reg},length(sig_ix),1);
            cpeaks = ctable.([creg 'peak'])(cdir_ix);

            if d == 1
                dir_lab = [roi_label{1} '_to_' roi_label{2}];
            elseif d == -1
                dir_lab = [roi_label{2} '_to_' roi_label{1}];
            else
                dir_lab = [roi_label{1} '_and_' roi_label{2}];
            end
 
            plot_name = [SBJ_id '_' model_id '_' an_id '_' stat_id '_' conn_id '_' rcn.plot_roi '_' ccat...
                '_' creg_name '_' rcn.hemi_str '_' rcn.view_str,'_' dir_lab];
            
            fig = fn_plot_recon_mesh(celec, roi_mesh, roi_mesh_lab, rcn, plot_name);
            %max_peak = max(abs(ctable.([creg 'peak'])));
            max_peak = mean(abs(ctable.([creg 'peak'])))*1.66;
            cclim = [-max_peak,max_peak];
            %cclim = [-max_peak-max_peak*0.01, max_peak+max_peak*0.01];
%             cclim = [mean(ctable.([creg 'peak'])(ctable.([creg 'peak'])<0))*1.25,...
%                      mean(ctable.([creg 'peak'])(ctable.([creg 'peak'])>0))*1.25];
%             if cmap_dir == -1
%                 ccmap = flipud(colormap(cmap));
%             else
%                 ccmap = colormap(cmap);
%             end
            ccmap = interp1([1 2 3], [cmap{2}; 255 255 255; cmap{1}]/255, linspace(1, 3, 200), 'linear');
            colvals = linspace(cclim(1),cclim(2),size(ccmap,1));
            %colormap(fig,ccmap)
            
            hold on;

            for pix = 1:size(cur_reg_elecs,1)
                if d == -1
                    dix1 = 2; dix2 = 1;
                else
                    dix1 = 1; dix2 = 2;
                end
                %ccolor =  colormap(reg_colors{reg};
                
                [~,ccolix] = min(abs(colvals-squeeze(cpeaks(pix))));
                ccolor = colors{cat};%ccmap(ccolix,:);%colors{cat};%
                pos1 = celec.chanpos(strcmp(celec.label,cur_reg_elecs{pix,dix1}),:);
                pos2 = celec.chanpos(strcmp(celec.label,cur_reg_elecs{pix,dix2}),:);
                l1 = plot3([pos1(1);pos2(1)],[pos1(2);pos2(2)],[pos1(3);pos2(3)],'-',...
                    'color',ccolor,'LineWidth',0.5);
                l1.Color(4) = 0.9;
                %         if dirs{reg}(pix) == 0
                %             sah = 0;
                %         else
                %             sah = 1;
                %         end
                %         [cx,cy,cz] = cylinder([1,0]);
                %         cz = cz*5;
                %         %xangle = acosd(max(min(dot([cx(1,end),cy(1,end)],pos2(1:2))/(norm([cx(1,end),cy(1,end)])*norm(pos2(1:2))),1),-1));
                %         %yangle = acosd(max(min(dot([cy(1,end),cz(1,end)],pos2(1:2))/(norm([cy(1,end),cz(1,end)])*norm(pos2(2:3))),1),-1));
                %         %M=makehgtform('translate',pos2*0.5);%,'xrotate',xangle,'yrotate',yangle);
                %         cx = cx + pos2(1)*0.5;
                %         cy = cy + pos2(2)*0.5;
                %         cz = cz + pos2(3)*0.5;
                %         surf(cx,cy,cz,'FaceColor',reg_colors{reg})%,'Parent',hgtransform('Matrix',M))
                % %         q = quiver3(pos1(1),pos1(2),pos1(3),pos2(1)-pos1(1),pos2(2)-pos1(2),pos2(3)-pos1(3),0,...
                % %             'color',reg_colors{reg},'ShowArrowHead',sah,'AutoScale','off','LineWidth',.75,...
                % %                 'MaxHeadSize',0.5);
                % %         disp(q)
                hold on
            end
            %colorbar
            ax = gca;
            %ax.CLim = cclim;
            if save_fig
                fig_fname = [out_dir, plot_name '.' fig_ftype];
                fig_fname = strrep(fig_fname,'*','x');
                if strcmp(fig_ftype,'fig')
                    saveas(fig,fig_fname);
                else
                    print(fig_fname,sprintf('-d%s',fig_ftype), sprintf('-r%d',fig_res));
                end
            end
            hold off
            clear ax
        end
    end
end
% %% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
% if save_fig
%     out_dir = [root_dir 'PRJ_Error/results/' an_dir '/GRP/recon_conn/' model_id '/' an_id '/' conn_id '/'];
%     if ~exist(out_dir,'dir'); [~] = mkdir(out_dir); end
% end
%
% for reg = 1:length(regs)
%     cdirs = unique(dirs{reg});
%     disp(cdirs)
%     for dd = 1:length(cdirs)
%         celec = elec;
%         d = cdirs(dd);
%         sig_ix = []; cdir_ix =[];
%         cdir_ix = find(ismember(dirs{reg},d));
%         cur_reg_elecs = sig_reg_elecs{reg}(cdir_ix,:);
%         sig_ix = find(ismember(celec.label, cur_reg_elecs));
%         celec.chanpos = celec.chanpos(sig_ix,:);
%         celec.chantype = celec.chantype(sig_ix,:);
%         celec.chanunit = celec.chanunit(sig_ix,:);
%         celec.label= celec.label(sig_ix,:);
%         celec.color = repmat(reg_colors{reg},length(sig_ix),1);
%         %celec.color = repmat([217 217 214]/255, size(celec.chanpos,1), 1);
%         %celec.color(sig_ix,:) = repmat(reg_colors{reg},length(sig_ix),1);
%
%         if d == 1
%             dir_lab = [roi_label{1} '_to_' roi_label{2}];
%         elseif d == -1
%             dir_lab = [roi_label{2} '_to_' roi_label{1}];
%         else
%             dir_lab = [roi_label{1} '_and_' roi_label{2}];
%         end
%         plot_name = [SBJ_id '_' model_id '_' an_id '_' conn_id '_' rcn.plot_roi '_' regs{reg}...
%             '_' rcn.hemi_str '_' rcn.view_str,'_' dir_lab];
%
%         fig = fn_plot_recon_mesh(celec, roi_mesh, roi_mesh_lab, rcn, plot_name);
%         max_peak = max(abs(peaks{reg}));
%         cclim = [-max_peak-max_peak*0.05, max_peak+max_peak*0.05];
%         ccmap = flipud(colormap(reg_cmaps{reg}));
%         colvals = linspace(cclim(1),cclim(2),size(ccmap,1));
%         colormap(fig,ccmap)
%
%         hold on;
%         for pix = 1:length(cur_reg_elecs)
%             if d == -1
%                 dix1 = 2; dix2 = 1;
%             else
%                 dix1 = 1; dix2 = 2;
%             end
%             %ccolor =  colormap(reg_colors{reg};
%             [~,ccolix] = min(abs(colvals-peaks{reg}(pix)));
%             ccolor = ccmap(ccolix,:);
%             pos1 = celec.chanpos(strcmp(celec.label,cur_reg_elecs{pix,dix1}),:);
%             pos2 = celec.chanpos(strcmp(celec.label,cur_reg_elecs{pix,dix2}),:);
%             l1 = plot3([pos1(1);pos2(1)],[pos1(2);pos2(2)],[pos1(3);pos2(3)],'-',...
%                 'color',ccolor,'LineWidth',0.5);
%             l1.Color(4) = 0.9;
%             %         if dirs{reg}(pix) == 0
%             %             sah = 0;
%             %         else
%             %             sah = 1;
%             %         end
%             %         [cx,cy,cz] = cylinder([1,0]);
%             %         cz = cz*5;
%             %         %xangle = acosd(max(min(dot([cx(1,end),cy(1,end)],pos2(1:2))/(norm([cx(1,end),cy(1,end)])*norm(pos2(1:2))),1),-1));
%             %         %yangle = acosd(max(min(dot([cy(1,end),cz(1,end)],pos2(1:2))/(norm([cy(1,end),cz(1,end)])*norm(pos2(2:3))),1),-1));
%             %         %M=makehgtform('translate',pos2*0.5);%,'xrotate',xangle,'yrotate',yangle);
%             %         cx = cx + pos2(1)*0.5;
%             %         cy = cy + pos2(2)*0.5;
%             %         cz = cz + pos2(3)*0.5;
%             %         surf(cx,cy,cz,'FaceColor',reg_colors{reg})%,'Parent',hgtransform('Matrix',M))
%             % %         q = quiver3(pos1(1),pos1(2),pos1(3),pos2(1)-pos1(1),pos2(2)-pos1(2),pos2(3)-pos1(3),0,...
%             % %             'color',reg_colors{reg},'ShowArrowHead',sah,'AutoScale','off','LineWidth',.75,...
%             % %                 'MaxHeadSize',0.5);
%             % %         disp(q)
%             hold on
%         end
%         colorbar
%         ax = gca;
%         ax.CLim = cclim;
%         if save_fig
%             fig_fname = [out_dir plot_name '.' fig_ftype];
%             fig_fname = strrep(fig_fname,'*','x');
%             if strcmp(fig_ftype,'fig')
%                 saveas(fig,fig_fname);
%             else
%                 print(fig_fname,sprintf('-d%s',fig_ftype), sprintf('-r%d',fig_res));
%             end
%         end
%         hold off
%         clear ax
%     end
% end
