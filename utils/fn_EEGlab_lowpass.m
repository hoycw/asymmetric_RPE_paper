function signal_lowpass = fn_EEGlab_lowpass(signal, srate, freq)
% Will use EEG lab to bandpass a signal
addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

signal_lowpass = eegfilt(signal, srate, [], freq);

end
