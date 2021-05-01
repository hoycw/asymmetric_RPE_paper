import numpy as np
import h5py
from scipy import signal
from scipy.stats.stats import pearsonr
from scipy.stats import zscore

import os


def columnwise_correlation(ypred, y, zscorea=True, zscoreb=True, axis=0):
    r'''Compute correlations efficiently
    Examples
    --------
    >>> x = np.random.randn(100,2)
    >>> y = np.random.randn(100,2)
    >>> cc = columnwise_correlation(x, y)
    >>> cc.shape
    (2,)
    >>> c1 = np.corrcoef(x[:,0], y[:,0])[0,1]
    >>> c2 = np.corrcoef(x[:,1], y[:,1])[0,1]
    >>> assert np.allclose(cc, np.r_[c1, c2])
    Notes
    -----
    Recall that the correlation cofficient is defined as
    .. math::
       \rho_{x, y} = \frac{cov(X,Y)}{var(x)var(y)}
    Since it is scale invariant, we can zscore and get the same
    .. math::
       \rho_{x, y} = \rho_{zscore(x), zscore(y)} = \frac{cov(X,Y)}{1*1} =
       \frac{1}{N}\frac{\sum_i^n \left(x_i - 0 \right) \left(y_i - 0 \right)}{1*1} =
       \frac{1}{N}\sum_i^n \left(x_i * y_i \right)
    '''
    if zscorea:
        y = zscore(y, axis=axis)
    if zscoreb:
        ypred = zscore(ypred, axis=axis)
    corr = (y * ypred).mean(axis=axis)
    return corr


def cpickle_dump(flname, datadict):
    '''
    Dump a dictionary containing some datas into a pickle file
    '''
    import cPickle
    assert isinstance(datadict, dict)
    print 'Saving %s...' % flname
    fl = open(flname, 'w')
    cPickle.dump(datadict, fl)
    fl.close()
    return


############### FOR NEURAL FEATURES #####################
not_hfb = range(0,7)

filepath = '/home/knight/cmerrick/notebooks_regression/mat_files/EC197/neural_data_mat.mat'
neural_mat = {}

with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        neural_mat[k] = np.array(v[not_hfb,:,:], dtype = np.float32)


neural_feat = neural_mat['neural_data_mat']

neural_feat = neural_feat.T

############# FOR KIN FEATURES ##############
kin_feat = np.hstack((1,4,5,6)) # movement + target model

# orignal
# kin lag mats
filepath = '/home/knight/cmerrick/notebooks_regression/mat_files/EC197/recon_lags_400.h5'
kin_lags_mat = {}
f = h5py.File(filepath)
with h5py.File(filepath,'r') as f:
    for k, v in f.items():
        kin_lags_mat[k] = np.array(v[kin_feat,:,:], dtype = np.float32)

kin_feat = kin_lags_mat['kinematics']
kin_feat = np.reshape(kin_feat,(kin_feat.shape[0]*kin_feat.shape[1],kin_feat.shape[2]))


kin_feat = kin_feat.T


# EC194 hand index
# cutting out long intervals between blocks
# ipsi_ind = np.hstack((range(5000,32000), range(37000,55000),range(58000,69685)))
# contra_ind = np.hstack((range(72800,101000), range(119800,147397)))

# # EC197
ipsi_ind = range(0,57117)
contra_ind = np.hstack(range(57117,114409))



# for testing ipsi and contra within hand ##

# lamdbas for ridge
l = np.logspace(0,8, num = 20)
num_elec = neural_feat.shape[1]
num_hands = 2
num_n_feat = neural_feat.shape[2]

# number of mutually exclusive test sets (each comprising of 20% of the data)
num_test_set = 5
test_set_run = [[0,1,2,3,4],[1,2,3,4,0],[2,3,4,0,1],[3,4,0,1,2],[4,0,1,2,3]]

# number of mutually exclusive validation sets within the estimation set (80% training, 20% validation)
num_val_fold = 5
val_fold_run = [[0,1,2,3,4],[1,2,3,4,0],[2,3,4,0,1],[3,4,0,1,2],[4,0,1,2,3]]

### correlation coefficient across folds and test sets ###
cc_test = np.zeros((num_hands,num_test_set,num_elec))

## to hold prediction test sets
prediction_test = {}
ecog_test = {}
feat_test = {}
feat_test_angle = {}
weights_reg = {}

n_feat = 4
# for n_feat in range(1):


for hand in range(num_hands):

    # init
    neural_mat = []
    kin_mat = []
    cc_val = np.zeros((num_test_set,num_val_fold,len(l),num_elec))

    if hand == 0:

        neural_mat = neural_feat[ipsi_ind,:,n_feat]
        kin_mat = kin_feat[ipsi_ind,:]

    elif hand == 1:

        neural_mat = neural_feat[contra_ind,:,n_feat]
        kin_mat = kin_feat[contra_ind,:]


    #### CREATING TEST SET ####
    for test_set in range(num_test_set):
    # for test_set in range(1):

        test_run = test_set_run[test_set]

        split_data_feat = np.array_split(kin_mat, num_test_set)
        split_data_ecog = np.array_split(neural_mat, num_test_set)

        ### for features/kinematics ###

        estimate_temp_f = []
        for ind in range(num_test_set-1):
            estimate_temp_f.extend(split_data_feat[test_run[ind]])

        test_temp_f = []
        test_temp_f.extend(split_data_feat[test_run[4]])

        estimate_feat = np.squeeze(estimate_temp_f)
        test_feat = np.squeeze(test_temp_f)

        ### for ecog ###

        estimate_temp_e = []
        for ind in range(num_test_set-1):
            estimate_temp_e.extend(split_data_ecog[test_run[ind]])

        test_temp_e = []
        test_temp_e.extend(split_data_ecog[test_run[4]])

        estimate_ecog = np.squeeze(estimate_temp_e)
        test_ecog = np.squeeze(test_temp_e)


        #### CREATING VALIDATION FOLD ####
        for fold in range(num_val_fold):

            val_fold = val_fold_run[fold]

            estimate_split_feat = np.array_split(estimate_feat, num_val_fold)
            estimate_split_ecog = np.array_split(estimate_ecog, num_val_fold)

            ### for features/kinematics ###

            train_temp_f = []
            for ind in range(num_val_fold-1):
                train_temp_f.extend(estimate_split_feat[val_fold[ind]])

            val_temp_f = []
            val_temp_f.extend(estimate_split_feat[val_fold[4]])

            train_feat = np.squeeze(train_temp_f)
            val_feat = np.squeeze(val_temp_f)

            ### for ecog ###

            train_temp_e = []
            for ind in range(num_val_fold-1):
                train_temp_e.extend(estimate_split_ecog[val_fold[ind]])

            val_temp_e = []
            val_temp_e.extend(estimate_split_ecog[val_fold[4]])

            train_ecog = np.squeeze(train_temp_e)
            val_ecog = np.squeeze(val_temp_e)


            # doing all computation we can outside lambda loop #
            XTX = np.dot(train_feat.T,train_feat)
            XTY = np.dot(train_feat.T,train_ecog)

            #### ridge regression ####

            # finding best lambda for each validation set
            for param in range(len(l)):
                XTXinv = np.linalg.inv(XTX +l[param]*np.identity(XTX.shape[0]))
                weights = np.dot(XTXinv,XTY)

                # create the prediction for the validation set
                # init
                num_timepnts = val_ecog.shape[0]
                num_elec = val_ecog.shape[1]
                prediction = np.zeros([num_timepnts])

                # record correlation from prediction set
                prediction = np.dot(val_feat,weights)
                mask = prediction[:,0] != 0
                cc_val[test_set,fold,param,:] = columnwise_correlation(prediction[mask,:],val_ecog[mask,:], zscorea=True, zscoreb=True, axis=0)


       # using best lambda to predict held-out test set
        for elec in range(test_ecog.shape[1]):

            max_corr = np.max(np.mean(cc_val[test_set,:,:,elec],0))

            best_lambda = np.squeeze(np.where(np.mean(cc_val[test_set,:,:,elec],0) == max_corr))

            XTXinv = np.linalg.inv(XTX + l[best_lambda]*np.identity(XTX.shape[0]))
            weights_reg[hand,test_set,elec] = np.dot(XTXinv,XTY[:,elec])

            # create the prediction for the test set
            prediction_test[hand,test_set,elec] = np.squeeze(np.sum([test_feat*weights_reg[hand,test_set,elec]],axis = 2))
            ecog_test[hand,test_set,elec] = test_ecog[:,elec]
            feat_test[hand,test_set,elec] = test_feat[:,200]
            feat_test_angle[hand,test_set,elec] = test_feat[:,600]

            # record correlation from prediction set
            mask = prediction_test[hand,test_set,elec] != 0
            cc_test[hand,test_set,elec] = columnwise_correlation(prediction_test[hand,test_set,elec][mask],ecog_test[hand,test_set,elec][mask], zscorea=True, zscoreb=True, axis=0)



## putting together all weights ##
num_elec = neural_feat.shape[1]
num_test_set = 5
num_weights = 1600
num_hand = 2
avg_weights = np.zeros((num_hand,num_test_set,num_elec,num_weights))


test_set_run = [1,2,3,4,0]

for hand in range(num_hand):

    for test_set in range(num_test_set):

        for elec in range(num_elec):

            avg_weights[hand,test_set,elec,:] = weights_reg[hand,test_set,elec]


## Make dictionary to save our vars ##
os.chdir('/home/knight/cmerrick/notebooks_regression/model_outputs_kin/EC197')

kin_dict = {}

kin_dict = {'avg_weights':avg_weights, 'cc_test':cc_test,'ecog_test':ecog_test, 'prediction_test':prediction_test}

kin_file_name = 'full_model_low_beta'
cpickle_dump(kin_file_name, kin_dict)

