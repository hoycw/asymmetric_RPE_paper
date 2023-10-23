# Asymmetric RPE paper code repository
This code repository is for all preprocessing, analysis, and plotting scripts related to Hoy, Quiroga-Martinez, et al. "Asymmetric coding of reward prediction errors in human insula and dorsomedial prefrontal cortex" (Nature Communications 2023)

## Data and Experimental Paradigm:
PsychoPy code for the paradigm is available at https://github.com/hoycw/target_time_scripts.git
All raw datasets are provided under Creative Commons Attrubition 4.0 license on Zenodo at https://doi.org/10.5281/zenodo.10023443.
This repo is citable via Zenodo: https://doi.org/10.5281/zenodo.10032478

## Manuscript: 
Manuscript has been accepted for publication at Nature Communications and is in press.
Initial biorxiv preprint: https://www.biorxiv.org/content/10.1101/2022.12.07.519496v1
Code written by Colin W. Hoy, David R. Quiroga-Martinez, and Eduardo Sandoval.

## Dependencies
OS: MacBook Pro running OS 10.13.6; MATLAB version R2017b; R (2022.12.0+353, lme4 package 1.1.31); Python 2.7 (not tested on any other platforms)
  - External toolboxes:
    - Fieldtrip (version d073bb2de): <http://www.fieldtriptoolbox.org/>
  - see colin_PRJ_Error_py2.7.yml for python environment dependencies, which can be replicated in a Conda (or similar) environment (see <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>)

## Overview of Code
Code execution to reproduce all analyses in the paper can be found in run_lists.
These run scripts contain parameters used to call each function. All individual analyses should run in under an hour, depending on computing capacity.
1. run01_preproc_pipeline.m
    - BHV scripts preprocess and analyze behavioral log files.
    - SBJ00, SBJ01, SBJ02, SBJ03, and SBJ04 scripts preprocess iEEG data, including analog event triggers.
2. run02_elec_pipeline.m
    - processes elec files from anatomical reconstruction pipeline to assign region-of-interest labels and plot on mesh surfaces
3. run03_signal_proc.m
    - SBJ03 scripts compute and plot difference waves.
4. run04_RL_LME.m
    - BHV02 scripts for RL model parameter estimation and plotting
    - SBJ08 scripts for HFA linear mixed models
5. run05_connectivity.m
    - SBJ10 and SBJ11 scripts for HFA connectivity modeling and plotting

## Execution Parameters
Specific parameters for different preprocessing and analysis scripts are loaded as options.
Generally, these options are written as executable MATLAB code that is specified when calling
and analysis script, which then run the relevant code to load the parameters inside that script.

- SBJ_lists: text lists of groups of subjects
  - 'preproc' is the cohort of 10 patients in this paper
- SBJ_vars: contains subject specific information
  - Raw data file names
  - iEEG information (e.g., channel labels)
  - Preprocessing information (e.g., bad channels, trials, etc.)
- proc_vars: preprocessing parameters for artifact rejection
- an_vars: analysis parameters for HFA computations
  - Channel selection, epoching, filtering, baseline correction, etc.
- stat_vars: statistics parameters for modeling HFA data
- plt_vars: plotting parameters
  - Epochs, styles, markers, significance, legends, etc.
- model_vars: RL parameters to model HFA in certain conditions
