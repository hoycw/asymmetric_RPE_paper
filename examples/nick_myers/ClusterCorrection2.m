function [p,praw] = ClusterCorrection2(dat, nSims, p_thresh)
%[p,praw] = ClusterCorrection2(dat, nSims, p_thresh)
%% Cluster Correct across any number of dimensions

if nargin < 3
    p_thresh = 0.05;
end
mincluster = 0;
[praw,tObs] = ttestfast(dat);
maxSize = zeros(nSims,1);
siz = size(tObs);
nsamp = size(dat,1);
indvec = {};
indvec(1:numel(siz)) = {':'};
for sim=1:nSims
    %permind = sign(rand(size(dat,1),1)-.5);
    %permind = repmat(permind,size(tObs));
    simDat  = dat;
    
    flipinds = rand(nsamp,1) < 0.5;
    indvec{1} = flipinds;
    simDat(indvec{:}) = -simDat(indvec{:});
    
    [p,tval] = ttestfast(simDat);
    
    CC = bwconncomp(p<=p_thresh);
    if CC.NumObjects == 0
        maxSize(sim) = 0;
    else        
        inds = CC.PixelIdxList;
        lens = cellfun(@numel, inds);
        inds = inds(lens > mincluster);
        if numel(inds) > 0
            tmpSize = zeros(numel(inds),1);
            for c=1:numel(inds)
                tmpSize(c) = abs(sum(tval(inds{c})));
            end    
            maxSize(sim) = max(tmpSize);
        else
            maxSize(sim) = 0;
        end
    end    
end

p = nan(size(praw));
CC = bwconncomp(praw<=p_thresh);
if CC.NumObjects ~= 0
    inds = CC.PixelIdxList;
    lens = cellfun(@numel, inds);
    inds = inds(lens > mincluster);        
    tmpSize = zeros(numel(inds),1);
    for c=1:numel(inds)
        ObsCluster = abs(sum(tObs(inds{c})));
        p(inds{c}) = (sum(maxSize>=ObsCluster)/nSims);
    end
end
end