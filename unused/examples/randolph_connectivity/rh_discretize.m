function bins = rh_discretize(X, Nbins, method)

% function to bin data for Mutual Information analysis
% use either UniCount (same N in every bin)
% or UniBins (linearly spaced bin sizes)

% equal counts method
if strcmp(method, 'UniCount') == 1

    [~,~,bins] = histcounts(X, [-inf, quantile(X, Nbins-1), inf]);

% requal bins methods
elseif strcmp(method, 'UniBins') == 1

    bins = discretize(X, linspace(min(X), max(X), Nbins+1));

else end

end