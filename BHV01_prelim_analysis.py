
# coding: utf-8

# In[16]:

from __future__ import division
# get_ipython().magic(u'matplotlib inline')
import sys 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as io
import pickle

import scipy.stats


# In[2]:

SBJ = sys.argv[1]#raw_input('Enter SBJ ID to process:')#'IR63'


# # Plotting Functions

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

def plot_tolerance_accuracy(SBJ, data, prdm, fig_fname):
    """Plot interval tolerance over experiment with accuracy overlay"""
    f, ax1 = plt.subplots()
    x = range(len(data))
    plot_title = '{0} Tolerance and Accuracy: easy={1:0.3f}; hard={2:0.3f}'.format(
                    SBJ, data[data['Condition']=='easy']['Hit'].mean(),
                    data[data['Condition']=='hard']['Hit'].mean())

    colors = {'easy': [0.5, 0.5, 0.5],#[c/255 for c in [77,175,74]],
              'hard': [1, 1, 1],#[c/255 for c in [228,26,28]],
              'accuracy': 'k'}#[c/255 for c in [55,126,184]]}
    scat_colors = {'easy': [1,1,1],#[c/255 for c in [77,175,74]],
              'hard': [0,0,0]}
    accuracy_colors = [scat_colors[accuracy.index[ix][1]] for ix in range(len(accuracy))]
    #scale = {'Hit Total': np.max(data['Tolerance'])/np.max(data['Hit Total']),
    #         'Score Total': np.max(data['Tolerance'])/np.max(data['Score Total'])}

    # Plot Tolerance Over Time
    ax1.plot(data['Tolerance'],'b',label='Tolerance')
    ax1.plot(x,[prdm['tol_lim'][0] for _ in x],'b--')
    ax1.plot(x,[prdm['tol_lim'][1] for _ in x],'b--')
    ax1.set_ylabel('Target Tolerance (s)', color='b')
    ax1.tick_params('y', colors='b')
    ax1.set_xlim([0,len(data)])
    ax1.set_ylim([0, 0.41])
    ax1.set_facecolor('white')
    ax1.grid(False)

    # Plot Accuracy per Block
    ax2 = ax1.twinx()
    # ax2.plot(data['Hit Total']/np.max(data['Hit Total']),'k',label='Hit Total')
    ax2.fill_between(x, 1, 0, where=data['Condition']=='easy',
                    facecolor=colors['easy'], alpha=0.3)#, label='hard')
    ax2.fill_between(x, 1, 0, where=data['Condition']=='hard',
                    facecolor=colors['hard'], alpha=0.3)#, label='easy')
    ax2.scatter(block_mid_ix, accuracy, s=50, c=accuracy_colors,
               edgecolors='k', linewidths=1)#colors['accuracy'])#,linewidths=2)
    ax2.set_ylabel('Accuracy', color=colors['accuracy'])
    ax2.tick_params('y', colors=colors['accuracy'])
    ax2.set_xlabel('Trials')
    ax2.set_xlim([0,len(data)])
    ax2.set_ylim([0, 1])
    ax2.set_facecolor('white')
    ax2.grid(False)

    plt.title(plot_title)

    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# In[6]:

def plot_ITI_histogram(SBJ, data, fig_fname):
    """Plot histogram of ITIs and save"""
    f,axes = plt.subplots(1,2)
    bins = np.arange(0,1.1,0.01)
    hist_real = sns.distplot(data['ITI'],bins=bins,kde=False,label=SBJ,ax=axes[0])
    hist_adj  = sns.distplot(data['ITI type'],bins=bins,kde=False,label=SBJ,ax=axes[1])
    axes[0].set_xlim([0, 1.1])
    axes[1].set_xlim([0, 1.1])
    plt.subplots_adjust(top=0.93)
    f.suptitle(SBJ)
    
    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# In[7]:

def plot_RT_histogram(SBJ, data, fig_fname):
    """Plot histogram of RTs and save"""
    f,ax = plt.subplots()
    hist = sns.distplot(data['RT'],label=SBJ)
    plt.subplots_adjust(top=0.9)
    hist.legend() # can also get the figure from plt.gcf()
    
    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# In[8]:

def plot_RT_by_ITI_hist(SBJ, data, fig_fname):
    """Plot RTs binned by ITI, with ANOVA stats for differences"""
    # Compute ANOVA for RT differences across ITI
    itis = np.unique(data['ITI type'])
    if len(prdm['ITIs'])==4:
        f,iti_p = scipy.stats.f_oneway(data.loc[data['ITI type']==itis[0],('RT')].values,
                                   data.loc[data['ITI type']==itis[1],('RT')].values,
                                   data.loc[data['ITI type']==itis[2],('RT')].values,
                                   data.loc[data['ITI type']==itis[3],('RT')].values)
    elif len(prdm['ITIs'])==3:
        f,iti_p = scipy.stats.f_oneway(data.loc[data['ITI type']==itis[0],('RT')].values,
                                   data.loc[data['ITI type']==itis[1],('RT')].values,
                                   data.loc[data['ITI type']==itis[2],('RT')].values)
    elif len(prdm['ITIs'])==2:
        f,iti_p = scipy.stats.ttest_ind(data.loc[data['ITI type']==itis[0],('RT')].values,
                                   data.loc[data['ITI type']==itis[1],('RT')].values)
    else:
        print 'WARNING: some weird paradigm version without 2, 3, or 4 ITIs!'
    
    # Plot RT histogram binned by ITI
    f, axes = plt.subplots(1,2)

    rt_bins = np.arange(0.7,1.3,0.01)
    for iti in itis:
        sns.distplot(data['RT'].loc[data['ITI type'] == iti],bins=rt_bins,label=str(round(iti,2)),ax=axes[0])
    axes[0].legend() # can also get the figure from plt.gcf()
    axes[0].set_xlim(min(rt_bins),max(rt_bins))

    # Factor Plot
    sns.boxplot(data=data,x='ITI type',y='RT',hue='ITI type',ax=axes[1])

    # Add overall title
    plt.subplots_adjust(top=0.9,wspace=0.3)
    f.suptitle(SBJ+' RT by ITI (p='+str(round(iti_p,4))+')') # can also get the figure from plt.gcf()

    # Save plot
    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# In[9]:

def plot_RT_adjustments(SBJ, data, fig_fname):
    """Plot RTs based on previous response timing (early/late) with t-test stats"""
    # t test for RT differences across ITI
    itis = np.unique(data['ITI type'])
    f,postlong_p = scipy.stats.ttest_ind(data.loc[data['postlong']==True,('dRT')].values,
                                data.loc[data['postlong']==False,('dRT')].values)
    
    # Plot RT histograms
    f, axes = plt.subplots(1,2)

    # RT Histogram
    drt_bins = np.arange(-0.6,0.6,0.025)
    sns.distplot(data['dRT'].loc[data['postlong']==True],bins=drt_bins,label='Post-Long',ax=axes[0])
    sns.distplot(data['dRT'].loc[data['postlong']==False],bins=drt_bins,label='Post-Short',ax=axes[0])
    axes[0].legend() # can also get the figure from plt.gcf()
    axes[0].set_xlim(min(drt_bins),max(drt_bins))

    # Factor Plot
    sns.boxplot(data=data,x='postlong',y='dRT',hue='postlong',ax=axes[1])

    # Add overall title
    plt.subplots_adjust(top=0.9,wspace=0.3)
    f.suptitle(SBJ+' RT by Previous RT Timing (p='+str(round(postlong_p,6))+')') # can also get the figure from plt.gcf()

    # Save plot
    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# In[10]:

def plot_RT_by_ITI_factor(SBJ, data, fig_fname):
    """Plot RT distribution factor plot by ITI type for easy and hard"""
    # RTs by condition
    # if len(prdm_params['ITIs'])==4:    # target_time v1.8.5+
    #     data['ITI type'] = ['short' if data['ITI'][ix]<0.5 else 'long' for ix in range(len(data))]
    #     ITI_plot_order = ['short','long']
    # elif len(prdm_params['ITIs'])==3:  # target_time v1.8.4 and below
    #     data['ITI type'] = ['short' if data['ITI'][ix]<prdm_params['ITI_bounds'][0] else 'long' \
    #                         if data['ITI'][ix]>prdm_params['ITI_bounds'][1] else 'medium'\
    #                         for ix in range(len(data))]
    #     ITI_plot_order = ['short','medium','long']
    # else:               # Errors for anything besides len(ITIs)==3,4
    #     assert len(prdm_params['ITIs'])==4

    plot = sns.factorplot(data=data,x='ITI type',y='dRT',hue='PE',col='Condition',kind='point',
                   ci=95);#,order=ITI_plot_order
    plt.subplots_adjust(top=0.9)
    plot.fig.suptitle(SBJ) # can also get the figure from plt.gcf()

    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# In[11]:

def plot_acc_by_ITI(SBJ, data, fig_fname):
    """Plot accuracy by ITI type for easy and hard"""

    # WARNING: I would need to go across subjects to get variance in accuracy by ITI
    plot = sns.factorplot(data=data,x='ITI type',y='Acc_ITI',col='Condition',kind='point',sharey=False,
                   ci=95);#,order=ITI_plot_order
    #plot.set(alpha=0.5)
    plt.subplots_adjust(top=0.9)
    plot.fig.suptitle(SBJ) # can also get the figure from plt.gcf()
    
    print 'Saving figure:', fig_fname
    plt.savefig(fig_fname)


# # Run and Save Plots

# In[13]:

# Set up directories
root_dir, app_dir = get_root_app_dir()

prj_dir     = os.path.join(root_dir,'PRJ_Error/')
results_dir = os.path.join(prj_dir,'results/')
data_dir    = os.path.join(prj_dir,'data/')
sbj_dir     = os.path.join(data_dir,SBJ)

fig_type = 'png'

# Get log names
logs = {}
with open(data_dir+'TT_behav_log_list.txt') as f:
    for line in f:
        (SBJ_id, log_names) = line.split(',',1)
        log_names = log_names.replace('\n','')
        logs[SBJ_id] = log_names.split(',')


# In[17]:

for log_ix in range(len(logs[SBJ])):
    # =====================================================================================================
    # LOAD DATA
    # =====================================================================================================
    if len(logs[SBJ])>1:
        run_suffix = '_R'+str(log_ix+1)
    else:
        run_suffix = ''

    # Load log
    log_filename = os.path.join(sbj_dir,'00_raw',logs[SBJ][log_ix])

    # Load Data
    behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav'+run_suffix+'.csv')
    print 'Processing log file for', SBJ, run_suffix, ':', behav_fname
    data = pd.read_csv(behav_fname)

    # Load Paradigm Variables
    prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars'+run_suffix+'.pkl')
    with open(prdm_fname, 'rb') as f:
        prdm = pickle.load(f)
    
    
    # =====================================================================================================
    # ADD ANALYSIS VARIABLES
    # =====================================================================================================
    # Label post-correct (PC), post-error (PE) trials
    data['PE'] = [False for _ in range(len(data))]
    for ix in range(len(data)):
        # Exclude training data and first trial of the block
        if (data.loc[ix,'Block']!=-1) and (data.loc[ix,'Trial']!=0):
            if data.loc[ix-1,'Hit']==0:
                data.loc[ix,'PE'] = True
    
    # Find middle of blocks to plot accuracy
    block_start_ix = data[data['Trial']==0].index
    block_mid_ix = [ix+prdm['n_trials']/2 for ix in block_start_ix[1:]]

    # Add in full_vis + E/H training: 0:4 + 5:19 = 10; 20:34 = 27.5 
    block_mid_ix.insert(0,np.mean([prdm['n_examples']+prdm['n_training'],
             prdm['n_examples']+2*prdm['n_training']]))   #examples
    block_mid_ix.insert(0,np.mean([0, prdm['n_examples']+prdm['n_training']]))
    #easy training (would be 12.5 if splitting examples/train)
    
    # Compute accuracy per block
    accuracy = data['Hit'].groupby([data['Block'],data['Condition']]).mean()
    acc_ITI = data['Hit'].groupby([data['ITI type'],data['Condition']]).mean()
    for ix in range(len(data)):
        data.loc[ix,'Accuracy'] = accuracy[data.loc[ix,'Block'],data.loc[ix,'Condition']]
        data.loc[ix,'Acc_ITI'] = acc_ITI[data.loc[ix,'ITI type'],data.loc[ix,'Condition']]
    
    # Break down by post-long and post-short trials
    data['postlong'] = [False if ix==0 else True if data['RT'].iloc[ix-1]>1 else False for ix in range(len(data))]

    # Compute change in RT
    data['dRT'] = [0 for ix in range(len(data))]
    for ix in range(len(data)-1):
        data.loc[ix+1,'dRT'] = data.loc[ix+1,'RT']-data.loc[ix,'RT']
    
    
    # =====================================================================================================
    # MAIN PLOT: Tolerance and Accuracy
    # =====================================================================================================
    fig_fname = os.path.join(results_dir,'BHV','tolerance',SBJ+'_tolerance'+run_suffix+'.'+fig_type)
    plot_tolerance_accuracy(SBJ, data, prdm, fig_fname)
    
    
    # =====================================================================================================
    # Remove examples and training data for remaining behavioral analyses
    # =====================================================================================================
    data_all = data
    # Exclude: Training/Examples, non-responses, first trial of each block
    if data[data['RT']<0].shape[0]>0:
        print 'WARNING: '+str(data[data['RT']<0].shape[0])+' trials with no response!'
    data = data[(data['Block']!=-1) & (data['RT']>0) & (data['ITI']>0)]
    
    
    # =====================================================================================================
    # ADDITIONAL PLOTS: RTs and Accuracy
    # =====================================================================================================
    # Plot ITI histogram
    fig_fname = os.path.join(results_dir,'BHV','ITIs',SBJ+'_ITI_hist'+run_suffix+'.'+fig_type)
    plot_ITI_histogram(SBJ, data, fig_fname)
    
    # Plot RT histogram
    fig_fname = os.path.join(results_dir,'BHV','RTs','histograms',SBJ+'_RT_hist'+run_suffix+'.'+fig_type)
    plot_RT_histogram(SBJ, data, fig_fname)
    
    # Plot RT histograms by ITIs
    fig_fname = os.path.join(results_dir,'BHV','RTs','hist_ITI',SBJ+'_RT_ITI_hist_box'+run_suffix+'.'+fig_type)
    plot_RT_by_ITI_hist(SBJ, data, fig_fname)
    
    # Plot RT adjustments (early/late)
    fig_fname = os.path.join(results_dir,'BHV','RTs','hist_dRT',SBJ+'_dRT_postlong_hist_box'+run_suffix+'.'+fig_type)
    plot_RT_adjustments(SBJ, data, fig_fname)
    
    # Plot RTs by ITI type (factor plot) within easy/hard
    fig_fname = os.path.join(results_dir,'BHV','RTs','hist_PE_ITI',SBJ+'_RT_PE_ITI_hit'+run_suffix+'.'+fig_type)
    plot_RT_by_ITI_factor(SBJ, data, fig_fname)
    
    # Plot Accuracy by ITI for easy/hard
    fig_fname = os.path.join(results_dir,'BHV','accuracy',SBJ+'_acc_ITI'+run_suffix+'.'+fig_type)
    plot_acc_by_ITI(SBJ, data, fig_fname)

    print


# ## Look for behavioral adjustments following short and long responses

# In[41]:

# plot = sns.factorplot(data=data_PL,x='ITI type',y='RT',hue='PE',col='Condition',kind='point',
#                ci=95,order=prdm['ITIs']);
# plt.subplots_adjust(top=0.9)
# plot.fig.suptitle(SBJ+'_post-long') # can also get the figure from plt.gcf()

# # plt.savefig(results_dir+'RT_plots/'+SBJ+'_RT_PE_ITI_hit'+fig_type)
# plot2 = sns.factorplot(data=data_PS,x='ITI type',y='RT',hue='PE',col='Condition',kind='point',
#                ci=95,order=prdm['ITIs']);
# plt.subplots_adjust(top=0.9)
# plot2.fig.suptitle(SBJ+'_post-short') # can also get the figure from plt.gcf()

# # plt.savefig(results_dir+'RT_plots/'+SBJ+'_RT_PE_ITI_hit'+fig_type)


# In[ ]:



