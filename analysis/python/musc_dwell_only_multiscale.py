'''
MuSC N-dimensional PFA
'''
import os
from pfa import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')
import numpy as np

celltype = 'musc'
exp1 = '20160626'
exp2 = '20160701_0'
exp3 = '20160701_1'
D = 2 # number of dimensions to consider
n_bins = 15 # number of bins for course-graining
log_trans = False # set to False for no log transformation
if log_trans:
    log_name = 'logTransformed'
else:
    log_name = ''

# Read in dataframes
def run(tau):
    fgf = []
    nofgf = []
    i = 1
    print('Importing data...')
    for exp in [exp1, exp2, exp3]:
        fgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'fgf2_exp_motility_statistics_split_' + str(tau).zfill(2) + '.csv'), add_plate=i))
        nofgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'nofgf2_exp_motility_statistics_split_' + str(tau).zfill(2) + '.csv'), add_plate=i))
        i += 1

    fgf_df = pd.concat(fgf)
    fgf_df.index = range(fgf_df.shape[0])
    nofgf_df = pd.concat(nofgf)
    nofgf_df.index = range(nofgf_df.shape[0])

    # Process dfs
    dx = [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5]
    binsx = [3, 5, 7, 10, 15, 2, 3, 5, 7, 10, 15, 20, 3, 5, 10, 3, 5]

    sig_array = np.zeros((len(dx), 6))
    for i in range(len(dx)):

        #fgf_stat = process_stationary(fgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_FGF2', D = dx[i], n_bins = binsx[i])
        #nofgf_stat = process_stationary(nofgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_noFGF2', D = dx[i], n_bins = binsx[i])

        # process dwell times
        plt.close('all')
        process_dwell_times(fgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_FGF2_'+str(dx[i])+'D'+'_'+str(binsx[i]).zfill(2)+'bins' + '_tau' + str(tau).zfill(2), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=15)
        plt.close('all')
        process_dwell_times(nofgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_noFGF2_'+str(dx[i])+'D'+'_'+str(binsx[i]).zfill(2)+'bins' + '_tau' + str(tau).zfill(2), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=15)
        plt.close('all')

for tau in [20, 25, 30]:
    run(tau)
