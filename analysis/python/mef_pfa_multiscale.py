'''
MEF N-dimensional PFA (multiscale)
'''
import os
from pfa import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')
import numpy as np

celltype = 'mef'
exp1wt = '20160711'
exp2wt = '20160925_0'
exp3wt = '20160925_1'
exp4wt = '20160927'

exp1mr = '20160917'
exp2mr = '20160918'

D = 2 # number of dimensions to consider
n_bins = 15 # number of bins for course-graining
log_trans = False # set to False for no log transformation
if log_trans:
    log_name = 'logTransformed'
else:
    log_name = ''

# Read in dataframes
def run(tau):
    wt = []
    i = 1
    print('Importing data...')
    for exp in [exp1wt, exp2wt, exp3wt, exp4wt]:
        wt.append(read_dataframe(os.path.join('../data/', celltype, 'wt', exp, 'exp_motility_statistics_split_' + str(tau) + '.csv'), add_plate=i))
        i += 1

    mycras = []
    i = 1
    for exp in [exp1mr, exp2mr]:
        mycras.append(read_dataframe(os.path.join('../data/', celltype, 'mycras', exp, 'exp_motility_statistics_split_' + str(tau) + '.csv'), add_plate=i))
        i += 1


    wt_df = pd.concat(wt)
    wt_df.index = range(wt_df.shape[0])
    mycras_df = pd.concat(mycras)
    mycras_df.index = range(mycras_df.shape[0])

    # Plot 2D scatter of 2D transition vectors
    wt_dv = process_dim_vectors(wt_df, os.path.join('../data/', celltype, 'nd_pfa', 'MEF_WT_' + str(D) + 'transition_vector_scatter_' + log_name + 'tau' + str(tau).zfill(2) + '.png'))
    mycras_dv = process_dim_vectors(mycras_df, os.path.join('../data/', celltype, 'nd_pfa', 'MEF_MycRas_' + str(D) + 'transition_vector_scatter_' + log_name + 'tau' + str(tau).zfill(2) + '.png'))

    # Save transition vectors for the 2D, (15,15) bin maps
    wt_sigp, wt_sv = process_samples(wt_df,  os.path.join('../data/', celltype, 'nd_pfa'), 'MEF_WT_plotSigT', D=2, n_bins=15, min_count=3, sv=True)
    wt_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MEF_WT_2DcgPFA_bins15_transition_vectors_' + 'tau' + str(tau).zfill(2) + '.csv'))
    mycras_sigp, mycras_sv = process_samples(mycras_df,  os.path.join('../data/', celltype, 'nd_pfa'), 'MEF_MycRas_plotSigT', D=2, n_bins=15, min_count=3, sv=True)
    mycras_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MEF_MycRas_2DcgPFA_bins15_transition_vectors_' + 'tau' + str(tau).zfill(2) + '.csv'))

    # Process dfs
    dx = [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5]
    binsx = [3, 5, 7, 10, 15, 2, 3, 5, 7, 10, 15, 20, 3, 5, 10, 3, 5]

    sig_array = np.zeros((len(dx), 6))
    for i in range(len(dx)):
        print('---')
        print('Processing WT...')
        print('D =',dx[i], '|', '# of bins = ', binsx[i])
        print('---')
        wt_sigp, wt_stat = process_samples(wt_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MEF_WT', D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=3)
        if np.any(wt_sigp):
            print('Plotting...')
            plt.close('all')
            f = plot_sig_transitions(wt_sigp,  os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2), 'MEF_WT_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '_tau' + str(tau).zfill(2) + '.png'), sigmarkers=True)
        print('---')
        print('Processing MycRas...')
        print('D =',dx[i], '|', '# of bins = ', binsx[i])
        print('---')
        mycras_sigp, mycras_stat = process_samples(mycras_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MEF_MycRas', D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=3)
        if np.any(mycras_sigp):
            print('Plotting...')
            plt.close('all')
            f = plot_sig_transitions(mycras_sigp,  os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2), 'MEF_MycRas_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name  + '_tau' + str(tau).zfill(2) + '.png'), sigmarkers=True)

        if np.any(wt_sigp):
            wtS = np.any(wt_sigp['reject_H0'])
        else:
            wtS = 0
        if np.any(mycras_sigp):
            mycrasS = np.any(mycras_sigp['reject_H0'])
        else:
            mycrasS = 0

        wtSs = wt_stat['p_val'].min() * 3 # Holm-Bonferroni correction
        mycrasSs = mycras_stat['p_val'].min() * 3 # Holm-Bonferroni correction

        sig_array[i,:] = [wtS, mycrasS, wtSs, mycrasSs, dx[i], binsx[i]]

        # process dwell times
        #plt.close('all')
        #process_dwell_times(fgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MEF_WT_'+str(dx[i])+'D'+'_'+str(binsx[i]).zfill(2)+'bins' + '_tau' + str(tau).zfill(2), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=20)
        #plt.close('all')
        #process_dwell_times(nofgf_df, os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2)), 'MEF_MycRas_'+str(dx[i])+'D'+'_'+str(binsx[i]).zfill(2)+'bins' + '_tau' + str(tau).zfill(2), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=20)

    sig_df = pd.DataFrame(sig_array)
    sig_df.columns = ['wt_sig', 'mycras_sig', 'wt_stat_p', 'mycras_stat_p', 'D', 'n_bins']
    sig_df.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'tau' + str(tau).zfill(2), 'MEF_Aggragate_Sig_ND_PFA' + log_name + '_tau' + str(tau).zfill(2) + '.csv'), sep=',', index=False)


for tau in [20, 25, 30]:
    run(tau)
