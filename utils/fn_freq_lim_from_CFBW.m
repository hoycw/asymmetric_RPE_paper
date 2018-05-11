function freq_lim = fn_freq_lim_from_CFBW(center_freq, bandwidth)
% Takes in a center frequency and bandwidth
% Returns a [2 1] array with the frequency range

freq_lim = [center_freq-(bandwidth/2) center_freq+(bandwidth/2)];

end
