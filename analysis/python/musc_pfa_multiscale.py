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

    # Plot 2D scatter of 2D transition vectors
    fgf_dv = process_dim_vectors(fgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau), 'MuSC_FGF2_' + str(D) + 'transition_vector_scatter_' + log_name + 'tau' + str(tau).zfill(2) + '.png'))
    nofgf_dv = process_dim_vectors(nofgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau), 'MuSC_noFGF2_' + str(D) + 'transition_vector_scatter_' + log_name + 'tau' + str(tau).zfill(2) + '.png'))

    # Save transition vectors for the 2D, (15,15) bin maps
    fgf_sigp, fgf_sv = process_samples(fgf_df,  os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau)), 'MuSC_FGF2_plotSigT', D=2, n_bins=15, min_count=3, sv=True)
    fgf_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MuSC_FGF2_2DcgPFA_bins15_transition_vectors_' + 'tau' + str(tau).zfill(2) + '.csv'))
    nofgf_sigp, nofgf_sv = process_samples(nofgf_df,  os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau)), 'MuSC_noFGF2_plotSigT', D=2, n_bins=15, min_count=3, sv=True)
    nofgf_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MuSC_noFGF2_2DcgPFA_bins15_transition_vectors_' + 'tau' + str(tau).zfill(2) + '.csv'))

    # Process dfs
    dx = [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5]
    binsx = [3, 5, 7, 10, 15, 2, 3, 5, 7, 10, 15, 20, 3, 5, 10, 3, 5]

    sig_array = np.zeros((len(dx), 6))
    for i in range(len(dx)):
        print('---')
        print('Processing FGF2+...')
        print('D =',dx[i], '|', '# of bins = ', binsx[i])
        print('---')
        fgf_sigp, fgf_stat = process_samples(fgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_FGF2', D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=3)
        print('Plotting...')
        plt.close('all')
        f = plot_sig_transitions(fgf_sigp,  os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2), 'MuSC_FGF2_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '_tau' + str(tau).zfill(2) + '.png'), sigmarkers=True)
        print('---')
        print('Processing FGF2-')
        print('D =',dx[i], '|', '# of bins = ', binsx[i])
        print('---')
        nofgf_sigp, nofgf_stat = process_samples(nofgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_noFGF2', D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=3)
        print('Plotting...')
        plt.close('all')
        f = plot_sig_transitions(nofgf_sigp,  os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2), 'MuSC_noFGF2_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name  + '_tau' + str(tau).zfill(2) + '.png'), sigmarkers=True)

        if np.any(fgf_sigp):
            fgfS = np.any(fgf_sigp['reject_H0'])
        else:
            fgfS = 0
        if np.any(nofgf_sigp):
            nofgfS = np.any(nofgf_sigp['reject_H0'])
        else:
            nofgfS = 0

        fgfSs = fgf_stat['p_val'].min() * 3 # Holm-Bonferroni correction
        nofgfSs = nofgf_stat['p_val'].min() * 3 # Holm-Bonferroni correction

        sig_array[i,:] = [np.any(fgf_sigp['reject_H0']), np.any(nofgf_sigp['reject_H0']), fgfSs, nofgfSs, dx[i], binsx[i]]

        # process dwell times
        plt.close('all')
        process_dwell_times(fgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_FGF2_'+str(dx[i])+'D'+'_'+str(binsx[i]).zfill(2)+'bins' + '_tau' + str(tau).zfill(2), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=20)
        plt.close('all')
        process_dwell_times(nofgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MuSC_noFGF2_'+str(dx[i])+'D'+'_'+str(binsx[i]).zfill(2)+'bins' + '_tau' + str(tau).zfill(2), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=20)

    sig_df = pd.DataFrame(sig_array)
    sig_df.columns = ['fgf_sig', 'nofgf_sig', 'fgf_stat_p', 'nofgf_stat_p', 'D', 'n_bins']
    sig_df.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2), 'MuSC_Aggragate_Sig_ND_PFA' + log_name + '_tau' + str(tau).zfill(2) + '.csv'), sep=',', index=False)


for tau in [20, 25, 30]:
    run(tau)
