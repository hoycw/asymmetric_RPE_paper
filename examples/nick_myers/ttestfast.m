function [p,tval] = ttestfast(x)
%function [p,tval] = ttestfast(x)

dim   = 1;
nans  = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end
df    = max(samplesize - 1,0);
%xmean = nanmean(x,dim);
%sdpop = nanstd(x,[],dim);
xmean = mean(x,dim);
sdpop = std(x,[],dim);
ser   = sdpop ./ sqrt(samplesize);
tval  = xmean ./ ser;
x     = tval.^2 ./ (df + tval.^2) ;
p     = 1 - betainc( x, 0.5, 0.5*df );
end