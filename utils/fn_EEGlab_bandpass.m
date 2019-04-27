function signal_bandpass = fn_EEGlab_bandpass(signal, srate, low_freq, hi_freq)
% Will use EEG lab to bandpass a signal
if exist('/home/knight/hoycw/','dir');app_dir='/home/knight/hoycw/Apps/';else app_dir='/Users/colinhoy/Code/Apps/';end
addpath(genpath([app_dir 'EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/']));

signal_highpass = eegfilt(signal, srate, low_freq, []);
signal_bandpass = eegfilt(signal_highpass, srate, [], hi_freq);


end
