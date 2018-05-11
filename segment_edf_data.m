clear
close all
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Edit Information in this section
data_dir = '/home/knight/hoycw/PRJ_Error/data/CP22/00_raw/'; %location where you want to save the experiment block
cd(data_dir);
path=strcat(data_dir,'ZA77717Y_1-1+.edf'); %location of the edf file 
my_block_name = 'CP22_raw.mat'; %task owner and task name

% start time of experimental block:
H_exp_start = 14;
M_exp_start = 06;
S_exp_start = 30;

% end time of experimental block
H_exp_end = 14;
M_exp_end = 25;
S_exp_end = 50;

%% Loading Data
%  Do NOT edit information in this section
disp('converting file from .edf to .mat...')

datafile = ft_read_data(path);
hdr = ft_read_header(path);
samplerate = hdr.Fs;
nSamples = hdr.nSamples
%% Segment Data
%  Do NOT edit information in this section
disp('segmenting file into specified block')
% start time of edf file 
H_block_start = hdr.orig.T0(4);
M_block_start = hdr.orig.T0(5);
S_block_start = hdr.orig.T0(6);

% end time of edf file 
H_block_end = floor(H_block_start + (nSamples/samplerate)/3600);
M_block_end = floor(M_block_start + ((nSamples/samplerate)-3600*(H_block_end-H_block_start))/60);
S_block_end = floor(S_block_start + ((nSamples/samplerate)-3600*(H_block_end-H_block_start)-60*(M_block_end-M_block_start)));
if S_block_end >= 60
    M_block_end = M_block_end + 1;
    S_block_end = S_block_end - 60;
end
if M_block_end >= 60
    H_block_end = H_block_end + 1;
    M_block_end = M_block_end - 60;
end

% double check that experiment start and end is 
% contained entirely in this file

if H_exp_end > H_block_end
    error('The experimental block extends beyond the data file.  Look for subsequent data files')
elseif H_exp_end == H_block_end && M_exp_end > M_block_end
    error('The experimental block extends beyond the data file.  Look for subsequent data files')
elseif H_exp_end == H_block_end && M_exp_end == M_block_end && S_exp_end > S_block_end
    error('The experimental block extends beyond the data file.  Look for subsequent data files')
end

if H_exp_start < H_block_start
    error('The experimental block starts before the data file')
elseif H_exp_start == H_block_start && M_exp_start < M_block_start
    error('The experimental block starts before the data file')
elseif H_exp_start == H_block_start && M_exp_start == M_block_start && S_exp_start < S_block_start
    error('The experimental block starts before the data file')
end


% how long into the block the experiment started 
H_start_diff = H_exp_start - H_block_start;
M_start_diff = M_exp_start - M_block_start;
if M_start_diff < 0
    M_start_diff = 60 + M_start_diff;
    H_start_diff = H_start_diff - 1;
end
S_start_diff = S_exp_start - S_block_start;
if S_start_diff < 0
    S_start_diff = 60 + S_start_diff;
    M_start_diff = M_start_diff - 1;
end

% how long into the block the experiment ended 
H_end_diff = H_exp_end - H_block_start;
M_end_diff = M_exp_end - M_block_start;
if M_end_diff < 0
    M_end_diff = 60 + M_end_diff;
    H_end_diff = H_end_diff - 1;
end
S_end_diff = S_exp_end - S_block_start;
if S_end_diff < 0
    S_end_diff = 60 + S_end_diff;
    M_end_diff = M_end_diff - 1;
end

% index data and create data file
n_exp_start = 1 + 3600*samplerate*H_start_diff + 60*samplerate*M_start_diff + samplerate*S_start_diff;
n_exp_end = 3600*samplerate*H_end_diff + 60*samplerate*M_end_diff + samplerate*S_end_diff;
data = datafile(:,n_exp_start:n_exp_end);

% edit header file
hdr.nSamples = n_exp_end - n_exp_start + 1;
header = hdr;

% save experiment
save(my_block_name, 'data', 'header')