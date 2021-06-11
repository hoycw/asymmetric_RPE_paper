function [sig_mask] = fn_threshold_pval(pvals,alpha,tail)
%% Apply one- or two-sided threshold to p values

if tail==1
    sig_mask = pvals<=alpha;
elseif tail==2
    sig_mask = pvals>=(1-alpha/2) | pvals<=alpha/2;
else
    error('Sided must be 1 or 2');
end

end