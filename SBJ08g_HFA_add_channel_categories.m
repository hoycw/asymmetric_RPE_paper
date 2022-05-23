function SBJ08g_HFA_add_channel_categories(an_id, model_id, stat_id, chpval_type)

if exist('/home/knight/hoycw/','dir'); root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
%% Load stats results:
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
fig_dir = [root_dir 'PRJ_Error/results/HFA/GRP/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
stats_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
load(stats_fname)

%% estimate categories
[cat_lab, ~, ~, ~, ~] = fn_puns_category_label_styles('puns');
for r = 1:numel(beta_chan.coefs)
    
    %array to store values
    beta_chan.chancat_ix{r} = cell(4,1);
    
    % Find pRPE and nRPE channels with any significant time points
    [pidx,~] = find(squeeze(beta_chan.([chpval_type 'vals']){r}(:,3,:) < .05));
    [nidx,~] = find(squeeze(beta_chan.([chpval_type 'vals']){r}(:,4,:) < .05));
    pidx = unique(pidx); nidx = unique(nidx);
    
    % find channels responding to either pRPE or nRPE (not both)
    pchan = pidx(~ismember(pidx,nidx));
    nchan = nidx(~ismember(nidx,pidx));
    
    % find channels responding to both pRPE and nRPE
    common_chan = pidx(ismember(pidx,nidx));
    
    % calculate channels responding to uRPE (slnchan) and sRPE (rwdchan)
    slnchan = [];
    rwdchan = [];
    for cch = 1:numel(common_chan)
        % select current common channel
        cchan = common_chan(cch);
        
        % get significant timepoints for pRPE and nRPE in this channel
        psig = squeeze(beta_chan.([chpval_type 'vals']){r}(cchan,3,:) < .05);
        nsig = squeeze(beta_chan.([chpval_type 'vals']){r}(cchan,4,:) < .05);
        
        % evaluate whether coefficients are +ve or -ve
        ppos = squeeze(double(beta_chan.coefs{r}(cchan,3,:) > 0));
        nneg = squeeze(double(beta_chan.coefs{r}(cchan,4,:) > 0));
        
        % if channel has only sig. positive coefs for pRPE and only sig.
        % negative coef for nRPE (or vice versa), then classify as
        % rwdchan (i.e. sRPE)
        if (~ismember(0,ppos(psig)) & ~ismember(1,nneg(nsig))) |...
                (~ismember(0,nneg(nsig)) & ~ismember(1,ppos(psig)))
            rwdchan = [rwdchan;cchan];
        else
            % if channel has any sig. negative coef for pRPE and any
            % sig. negative coef for nRPE, then clasify as sRPE ONLY
            % if not at same time points.
            
            %common significant times
            %commtidx = find(double(psig).*double(nsig));
            
            %compare coefficients
            %if sum(ppos(commtidx) ~= nneg(commtidx)) > 0
                %rwdchan = [rwdchan;cchan];
            %else
                % Everyhting else is uRPE
            slnchan = [slnchan;cchan];
            %end
        end
    end
    
    cat_struc = [];
    cat_struc.pRPE = pchan;
    cat_struc.nRPE = nchan;
    cat_struc.sRPE = rwdchan;
    cat_struc.uRPE = slnchan;
    
    for ctg = 1:numel(cat_lab)
        beta_chan.chancat_ix{r}{ctg} = cat_struc.(cat_lab{ctg});
    end
end
beta_chan.chancat_label = cat_lab;
beta_chan.pval_type = chpval_type;
save(stats_fname,'-v7.3', 'beta_chan')