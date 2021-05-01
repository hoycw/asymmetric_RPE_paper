
# coding: utf-8

# In[1]:

# get_ipython().magic(u'matplotlib inline')
import sys 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as io
import pickle
import json


# In[2]:

SBJ = sys.argv[1]#raw_input('Enter SBJ ID to process:')#'IR63'


# # Log Processing Functions

# In[3]:

def get_root_app_dir():
    """Identify the correct root and app directories for HWNI cluster or local"""
    if os.path.isdir('/home/knight/'):
        root_dir = '/home/knight/hoycw/'
        app_dir  = root_dir + 'Apps/'
    elif os.path.isdir('/Volumes/hoycw_clust/'):
        root_dir = '/Volumes/hoycw_clust/'
        app_dir  = '/Users/colinhoy/Code/Apps/'

    # Raise NameError if root_dir was not assigned
    root_dir
    
    return root_dir, app_dir


# In[4]:

def parse_prdm_vars(log):
    """Parse Target Time log lines to obtain paradigm design varaiables"""
    prdm = {}
    for line in log:
        # Script version
        if line.find('paradigm_name =')!=-1:
            prdm['prdm_name'] = line[line.find('= ')+2:line.find('\n')]
        if line.find('paradigm_version =')!=-1:
            prdm['prdm_version'] = line[line.find('= ')+2:line.find('\n')]

        # Timing variables
        if line.find('interval_dur =')!=-1:
            prdm['target'] = float(line[line.find('= ')+2:line.find('\n')])
        if line.find('feedback_delay =')!=-1:
            prdm['fb_delay'] = float(line[line.find('= ')+2:line.find('\n')])
        if line.find('feedback_dur =')!=-1:
            prdm['fb'] = float(line[line.find('= ')+2:line.find('\n')])
        if line.find('trigger_dur =')!=-1:
            prdm['trig_dur'] = float(line[line.find('= ')+2:line.find('\n')])

        # ITIs and boundaries between them
        if line.find('ITIs')!=-1:
            prdm['ITIs'] = [float(string)                         for string in line[line.find('[')+1:line.find(']')].split(',')]

        # Tolerance limits/clamps
        if line.find('tolerance_lim')!=-1:
            prdm['tol_lim'] = [float(string)                              for string in line[line.find('[')+1:line.find(']')].split(',')]

        # Trial count variables
        if line.find('n_blocks')!=-1:
            prdm['n_blocks'] = int(line[line.find('=')+2:])
        if line.find('n_trials')!=-1:
            prdm['n_trials'] = int(line[line.find('=')+2:])
        if line.find('n_examples')!=-1:
            prdm['n_examples'] = int(line[line.find('=')+2:])
        elif line.find('n_fullvis')!=-1:
            prdm['n_examples'] = int(line[line.find('=')+2:])
        if line.find('n_training')!=-1:
            prdm['n_training'] = int(line[line.find('=')+2:])

    # Add missing items from early log files
    if 'prdm_name' not in prdm:
        prdm['prdm_name'] = 'not_logged'
    if 'prdm_version' not in prdm:
        prdm['prdm_version'] = '<1.8.5'
    if 'n_examples' not in prdm:
        prdm['n_examples'] = int(-1)
    if 'n_training' not in prdm:
        prdm['n_training'] = int(-1)

    prdm['trl_len'] = prdm['target']+                prdm['fb_delay']+prdm['fb']
    
    return prdm


# In[19]:

def extract_trl_info(log,prdm):
    """Extract trial info from Target Time behavioral log file"""
    # Separate informative lines
    resp_lines = [line for line in log if line.find('Outcome=')!=-1]
    
    # Extract info from lines
    data = pd.DataFrame({'Block': [line[line.find('B')+1] for line in resp_lines],
                         'Trial': [int(line[line.find('_T')+2:line.find(':')]) for line in resp_lines],
                         'Feedback': ['W' if line.count('WIN')>0 else \
                                      'L' if line.count('LOSE')>0 else \
                                      'N' for line in resp_lines],
                         'RT': [line[line.find('RT')+5:line.find('RT')+5+13].strip() for line in resp_lines],
                         'Tolerance': [float(line[line.find('tol')+12:line.find('\n')]) for line in resp_lines],
                         'Timestamp': [float(line[:line.find('.')+4]) for line in resp_lines]
                        })
    
    # Add Sounds
    if prdm['prdm_version'][0]=='1':
        # Only win or loss in v1.*
        data['Sound'] = ['new_win_sound' if data['Feedback'][ix]=='W' else 'new_loss_sound' for ix in range(len(data))]
    else:
        # win, loss, or oddball in v2.*
        sound_lines = [line for line in log if line.find('SOUND =')!=-1]
        data['Sound'] = [line[line.find('sounds/')+7:line.find('.wav')] for line in sound_lines]
    
    # Fix Reversals, Block, and missed RTs
    for ix in range(len(data)):
        # Fix Surprise sounds
        if data['Sound'][ix].find('\\')!=-1:
            data.loc[ix,'Sound'] = data['Sound'][ix][data['Sound'][ix].find('\\')+1:]
        # Fix RTs
        if data['RT'][ix][0:2]=='-1':#No response
            data.loc[ix,'RT'] = -1
            #data.loc[ix,'Score'] = 0            # !!! may change depending on version !!!!
        else:# Real Responses
            if data['RT'][ix].find(';')!=-1:# shorter number of digits, clip ';'
                data.loc[ix,'RT'] = float(data['RT'][ix][:data['RT'][ix].find(';')])
            else:
                data.loc[ix,'RT'] = float(data['RT'][ix])

        # Fix Block coding
        if data['Block'][ix]=='T':# Training
            data.loc[ix,'Block'] = -1
        else:
            data.loc[ix,'Block'] = int(data['Block'][ix])

    data['Hit'] = [1 if abs(data['RT'][ix]-prdm['target']) <= data['Tolerance'][ix] else 0 for ix in range(len(data))]
    data['Score'] = [0 if data['Feedback'][ix]=='N' or data['RT'][ix]== -1 else                      100 if data['Feedback'][ix]=='W' else -100 for ix in range(len(data))]

    # Mark trials with bad feedback (apparently my logic has error up to 10 ms... damn it!)
    data['bad_fb'] = [False if data['Feedback'][ix]=='N' or data['RT'][ix]==-1 or                       (data['Feedback'][ix]=='W' and data['Hit'][ix]==1) or                       (data['Feedback'][ix]=='L' and data['Hit'][ix]==0)                       else True for ix in range(len(data))]

    # Add condition based on logging version
    if prdm['prdm_version'][0]=='2':
        data['Condition'] = [line[line.find('condition')+12:line.find('condition')+16] for line in resp_lines]
    else:
        data['Condition'] = [line[line.find('_type')+8:line.find('_type')+12] for line in resp_lines]

    # Calculate ITIs
    ITI_bounds = np.mean([prdm['ITIs'][:-1], prdm['ITIs'][1:]],0)
    data['ITI'] = [data['Timestamp'][ix]-data['Timestamp'][ix-1]-prdm['trl_len'] if ix!=0 else 0                    for ix in range(len(data))]
    data.loc[data['Trial']==0,'ITI'] = 0
    data.loc[data['Block']==-1,'ITI'] = 0
    # Match real ITIs to ITI categories
    ITI_bin_edges = np.insert(ITI_bounds, ITI_bounds.shape[0], prdm['ITIs'][-1]+0.1)
    ITI_bin_edges = np.insert(ITI_bin_edges, 0, 0)
    data['ITI type'] = [prdm['ITIs'][np.argmax(np.histogram(data['ITI'][ix],bins=ITI_bin_edges)[0])]            if data['ITI'][ix]!=0 else 0 for ix in range(len(data))]
    # if len(prdm['ITIs'])==4:    # target_time v1.8.5+
    # elif len(prdm['ITIs'])==3:  # target_time v1.8.4 and below
    # else:               # Errors for anything besides len(ITIs)==3,4
    #     assert len(prdm['ITIs'])==4

    # Print stats on bad feedback
    if any(data['bad_fb']):
        bad_ix = [ix for ix in range(len(data)) if data['bad_fb'][ix]]
        tmp = data.ix[bad_ix]
        tmp['error'] = abs(tmp['RT']-prdm['target'])
        tmp['error-tol'] = tmp['Tolerance']-tmp['error']
        print 'WARNING!!! Bad feedback found on {0} trials!'.format(len(bad_ix))
        print 'Max diff between RT error and tolerance = {0}, mean = {1}'.format(                                                    np.max(tmp['error-tol']),np.mean(np.abs(tmp['error-tol'])))
        tmp.ix[:,{'Tolerance','error','WIN','Hit','err-tol'}]   
    
    return data


# # Process Logs

# In[6]:

# Set up directories
root_dir, app_dir = get_root_app_dir()
prj_dir = os.path.join(root_dir,'PRJ_Error/')
results_dir = os.path.join(prj_dir,'results/')
fig_type = '.png'
data_dir = os.path.join(prj_dir,'data/')
sbj_dir  = os.path.join(data_dir,SBJ)

# Get log names
logs = {}
with open(data_dir+'TT_behav_log_list.txt') as f:
    for line in f:
        (SBJ_id, log_names) = line.split(',',1)
        log_names = log_names.replace('\n','')
        logs[SBJ_id] = log_names.split(',')


# In[20]:

# Process Logs
for log_ix in range(len(logs[SBJ])):
    if len(logs[SBJ])>1:
        run_suffix = '_R'+str(log_ix+1)
    else:
        run_suffix = ''

    # Load log
    log_filename = os.path.join(sbj_dir,'00_raw',logs[SBJ][log_ix])
    print 'Processing log file for', SBJ, run_suffix, ':', log_filename

    log_file = open(log_filename,'r')
    log = log_file.readlines()
    log_file.close()

    # Parse paradigm variables
    prdm = parse_prdm_vars(log)

    # Save paradigm variables
    prdm_fname = os.path.join(root_dir,'PRJ_Error','data',SBJ,'03_events',SBJ+'_prdm_vars'+run_suffix+'.pkl')
    print 'Saving paradigm variables file: ', prdm_fname
    with open(prdm_fname, 'wb') as f: # Python readable
        pickle.dump(prdm, f, pickle.HIGHEST_PROTOCOL)

    prdm_fname = prdm_fname.replace('.pkl','.mat')
    print 'Saving paradigm variables file: ', prdm_fname
    io.savemat(prdm_fname,prdm) # MATLAB Readable

    # Parse trial info
    data = extract_trl_info(log,prdm)

    # Save trial info
    behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav'+run_suffix+'.csv')
    print 'Saving trial info file: ', behav_fname
    data.to_csv(behav_fname,index_label='Total_Trial')
    print


# In[ ]:



