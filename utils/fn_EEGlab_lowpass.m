function signal_lowpass = fn_EEGlab_lowpass(signal, srate, freq)
% Will use EEG lab to bandpass a signal
if exist('/home/knight/hoycw/','dir');app_dir='/home/knight/hoycw/Apps/';else app_dir='/Users/colinhoy/Code/Apps/';end
addpath(genpath([app_dir 'EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/']));

signal_lowpass = eegfilt(signal, srate, [], freq);

end
