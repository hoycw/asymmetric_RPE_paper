function [eeg_data, Sorig, Sclean, f, amps, freqs, g] = cleanline2(eeg_data, s_rate, varargin)
% CWH: This modifided version from KLA hopefully runs without EEGLab and taking just [trials,time] matrix
% Mandatory             Information
% --------------------------------------------------------------------------------------------------
% EEG                   EEGLAB data structure
% --------------------------------------------------------------------------------------------------
%
% Optional              Information
% --------------------------------------------------------------------------------------------------
% LineFrequencies:      Line noise frequencies to remove                                                                      
%                       Input Range  : Unrestricted                                                                           
%                       Default value: 60  120                                                                                
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% ScanForLines:         Scan for line noise                                                                                   
%                       This will scan for the exact line frequency in a narrow range around the specified LineFrequencies    
%                       Input Range  : Unrestricted                                                                           
%                       Default value: 1                                                                                      
%                       Input Data Type: boolean                                                                              
%                                                                                                                             
% LineAlpha:            p-value for detection of significant sinusoid                                                                        
%                       Input Range  : [0  1]                                                                                 
%                       Default value: 0.01                                                                                   
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% Bandwidth:            Bandwidth (Hz)                                                                                        
%                       This is the width of a spectral peak for a sinusoid at fixed frequency. As such, this defines the     
%                       multi-taper frequency resolution.                                                                     
%                       Input Range  : Unrestricted                                                                           
%                       Default value: 1                                                                                      
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% SignalType:          Type of signal to clean                                                                               
%                       Cleaned ICA components will be backprojected to channels. If channels are cleaned, ICA activations    
%                       are reconstructed based on clean channels.                                                            
%                       Possible values: 'Components','Channels'                                                              
%                       Default value  : 'Components'                                                                         
%                       Input Data Type: string                                                                               
%                                                                                                                             
% ChanCompIndices:      IDs of Chans/Comps to clean                                                                           
%                       Input Range  : Unrestricted                                                                           
%                       Default value: 1:152                                                                                  
%                       Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                             
% SlidingWinLength:     Sliding window length (sec)                                                                           
%                       Default is the epoch length.                                                                          
%                       Input Range  : [0  4]                                                                                 
%                       Default value: 4                                                                                      
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% SlidingWinStep:       Sliding window step size (sec)                                                                        
%                       This determines the amount of overlap between sliding windows. Default is window length (no           
%                       overlap).                                                                                             
%                       Input Range  : [0  4]                                                                                 
%                       Default value: 4                                                                                      
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% SmoothingFactor:      Window overlap smoothing factor                                                                       
%                       A value of 1 means (nearly) linear smoothing between adjacent sliding windows. A value of Inf means   
%                       no smoothing. Intermediate values produce sigmoidal smoothing between adjacent windows.               
%                       Input Range  : [1  Inf]                                                                               
%                       Default value: 100                                                                                    
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% PaddingFactor:        FFT padding factor                                                                                    
%                       Signal will be zero-padded to the desired power of two greater than the sliding window length. The    
%                       formula is NFFT = 2^nextpow2(SlidingWinLen*(PadFactor+1)). e.g. For SlidingWinLen = 500, if PadFactor = -1, we    
%                       do not pad; if PadFactor = 0, we pad the FFT to 512 points, if PadFactor=1, we pad to 1024 points etc.                                                                                                  
%                       Input Range  : [-1  Inf]                                                                              
%                       Default value: 2                                                                                      
%                       Input Data Type: real number (double)                                                                 
%                                                                                                                             
% ComputeSpectralPower: Visualize Original and Cleaned Spectra                                                                
%                       Original and clean spectral power will be computed and visualized at end                              
%                       Input Range  : Unrestricted                                                                           
%                       Default value: true                                                                                      
%                       Input Data Type: boolean                                                                              
%                                                                                                                             
% NormalizeSpectrum:    Normalize log spectrum by detrending (not generally recommended)                                                                     
%                       Input Range  : Unrestricted                                                                           
%                       Default value: 0                                                                                      
%                       Input Data Type: boolean                                                                              
%                                                                                                                             
% VerboseOutput:        Produce verbose output                                                                                
%                       Input Range  : [true false]                                                                           
%                       Default value: true                                                                                      
%                       Input Data Type: boolean                                                                
%                                                                                                                             
% PlotFigures:          Plot Individual Figures                                                                               
%                       This will generate figures of F-statistic, spectrum, etc for each channel/comp while processing       
%                       Input Range  : Unrestricted                                                                           
%                       Default value: 0                                                                                      
%                       Input Data Type: boolean  
%
% --------------------------------------------------------------------------------------------------
% Output                Information
% --------------------------------------------------------------------------------------------------
% EEG                   Cleaned EEG dataset
% Sorig                 Original multitaper spectrum for each component/channel
% Sclean                Cleaned multitaper spectrum for each component/channel
% f                     Frequencies at which spectrum is estimated in Sorig, Sclean
% amps                  Complex amplitudes of sinusoidal lines for each
%                       window (line time-series for window i can be
%                       reconstructed by creating a sinudoid with frequency f{i} and complex 
%                       amplitude amps{i})
% freqs                 Exact frequencies at which lines were removed for
%                       each window (cell array)
% g                     Parameter structure. Function call can be
%                       replicated exactly by calling >> cleanline(EEG,g);
%
% Usage Example:
% EEG = pop_cleanline(EEG, 'Bandwidth',2,'ChanCompIndices',[1:EEG.nbchan],                  ...
%                          'SignalType','Channels','ComputeSpectralPower',true,             ...
%                          'LineFrequencies',[60 120] ,'NormalizeSpectrum',false,           ...
%                          'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,          ...
%                          'ScanForLines',true,'SmoothingFactor',100,'VerboseOutput',1,    ...
%                          'SlidingWinLength',EEG.pnts/EEG.srate,'SlidingWinStep',EEG.pnts/EEG.srate);
%
% See Also:
% pop_cleanline()

% Author: Tim Mullen, SCCN/INC/UCSD Copyright (C) 2011
% Date:   Nov 20, 2011
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

EEG.nbchan = size(eeg_data,1);
EEG.pnts = size(eeg_data,2);
EEG.trials = 1;
EEG.srate = s_rate;

g = arg_define([0 1], varargin, ...
    arg({'linefreqs','LineFrequencies'},[60 120],[],'Line noise frequencies to remove.'),...
    arg({'scanforlines','ScanForLines'},true,[],'Scan for line noise. This will scan for the exact line frequency in a narrow range around the specified LineFrequencies'),...
    arg({'p','LineAlpha','alpha'},0.01,[0 1],'p-value for detection of significant sinusoid'), ...
    arg({'bandwidth','Bandwidth'},2,[],'Bandwidth (Hz). This is the width of a spectral peak for a sinusoid at fixed frequency. As such, this defines the multi-taper frequency resolution.'), ...
    arg({'chanlist','ChanCompIndices','ChanComps'},sprintf('1:%d',EEG.nbchan),[1 EEG.nbchan],'Indices of Channels/Components to clean.','type','expression'),...
    arg({'winsize','SlidingWinLength'},fastif(EEG.trials==1,4,EEG.pnts/EEG.srate),[0 EEG.pnts/EEG.srate],'Sliding window length (sec). Default for epoched data is the epoch length. Default for continuous data is 4 seconds'), ...
    arg({'winstep','SlidingWinStep'},fastif(EEG.trials==1,1,EEG.pnts/EEG.srate),[0 EEG.pnts/EEG.srate],'Sliding window step size (sec). This determines the amount of overlap between sliding windows. Default for epoched data is window length (no overlap). Default for continuous data is 1 second.'), ...
    arg({'tau','SmoothingFactor'},100,[1 Inf],'Window overlap smoothing factor. A value of 1 means (nearly) linear smoothing between adjacent sliding windows. A value of Inf means no smoothing. Intermediate values produce sigmoidal smoothing between adjacent windows.'), ...
    arg({'pad','PaddingFactor'},2,[-1 Inf],'FFT padding factor. Signal will be zero-padded to the desired power of two greater than the sliding window length. The formula is NFFT = 2^nextpow2(SlidingWinLen*(PadFactor+1)). e.g. For N = 500, if PadFactor = -1, we do not pad; if PadFactor = 0, we pad the FFT to 512 points, if PadFactor=1, we pad to 1024 points etc.'), ...
    arg({'computepower','ComputeSpectralPower'},true,[],'Visualize Original and Cleaned Spectra. Original and clean spectral power will be computed and visualized at end'),  ...
    arg({'normSpectrum','NormalizeSpectrum'},false,[],'Normalize log spectrum by detrending. Not generally recommended.'), ...
    arg({'verb','VerboseOutput','VerbosityLevel'},true,[],'Produce verbose output.'), ...
    arg({'plotfigures','PlotFigures'},false,[],'Plot Individual Figures. This will generate figures of F-statistic, spectrum, etc for each channel/comp while processing') ...
    );
      
% defaults
[Sorig, Sclean, f, amps, freqs] = deal([]);

% set up multi-taper parameters
hbw             = g.bandwidth/2;   % half-bandwidth
params.tapers   = [hbw, g.winsize, 1];
params.Fs       = s_rate;
params.g.pad    = g.pad;
movingwin       = [g.winsize g.winstep];

% NOTE: params.tapers = [W, T, p] where:
% T==frequency range in Hz over which the spectrum is maximally concentrated 
%    on either side of a center frequency (half of the spectral bandwidth)
% W==time resolution (seconds)
% p is used for num_tapers = 2TW-p (usually p=1).

SlidingWinLen = movingwin(1)*params.Fs;
if params.g.pad>=0
    NFFT = 2^nextpow2(SlidingWinLen*(params.g.pad+1));
else
    NFFT = SlidingWinLen;
end

k=0;
for ch=g.chanlist
    
    % extract data as [chans x frames*trials]
    data = squeeze(eeg_data(ch,:));
    
    % DO THE MAGIC!
    [datac,datafit,amps,freqs]=rmlinesmovingwinc(data,movingwin,g.tau,params,g.p,fastif(g.plotfigures,'y','n'),g.linefreqs,fastif(g.scanforlines,params.tapers(1),[]));   
    
    % append to clean dataset any remaining samples that were not cleaned 
    % due to sliding window and step size not dividing the data length
    ndiff = length(data)-length(datac);
    if ndiff>0
        datac(end:end+ndiff) = data(end-ndiff:end);
    end
        
    eeg_data(ch,:) = datac';

end


