'''
MEF N-dimensional PFA
'''
import os
from pfa import *
import pandas as pd
import numpy as np

celltype = 'mef'
exp1wt = '20160711'
exp2wt = '20160925_0'
exp3wt = '20160925_1'
exp4wt = '20160927'

exp1mr = '20160917'
exp2mr = '20160918'

D = 2 # number of dimensions to consider
n_bins = 3 # number of bins for course-graining
log_trans = False # set to None for no log transformation
if log_trans:
    log_name = 'logTransformed'
else:
    log_name = ''

# Read in dataframes

wt = []
i = 1
print('Importing data...')
for exp in [exp1wt, exp2wt, exp3wt, exp4wt]:
    wt.append(read_dataframe(os.path.join('../data/', celltype, 'wt', exp, 'exp_motility_statistics_split_20.csv'), add_plate=i))
    i += 1

mycras = []
i = 1
for exp in [exp1mr, exp2mr]:
    mycras.append(read_dataframe(os.path.join('../data/', celltype, 'mycras', exp, 'exp_motility_statistics_split_20.csv'), add_plate=i))
    i += 1


wt_df = pd.concat(wt)
wt_df.index = range(wt_df.shape[0])
mycras_df = pd.concat(mycras)
mycras_df.index = range(mycras_df.shape[0])

# Plot 2D scatter of 2D transition vectors
wt_dv = process_dim_vectors(wt_df, os.path.join('../data/', celltype, 'nd_pfa', 'MEF_WT_' + str(D) + 'transition_vector_scatter' + log_name + '.png'))
mycras_dv = process_dim_vectors(mycras_df, os.path.join('../data/', celltype, 'nd_pfa', 'MEF_MycRas_' + str(D) + 'transition_vector_scatter' + log_name + '.png'))

# Save transition vectors for the 2D, (15,15) bin maps
wt_sigp, wt_sv = process_samples(wt_df,  os.path.join('../data/', celltype, 'nd_pfa'), 'MEF_WT_plotSigT', D=2, n_bins=15, min_count=3, sv=True)
wt_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MEF_WT_2DcgPFA_bins15_transition_vectors.csv'))
mycras_sigp, mycras_sv = process_samples(mycras_df,  os.path.join('../data/', celltype, 'nd_pfa'), 'MEF_MycRas_plotSigT', D=2, n_bins=15, min_count=3, sv=True)
mycras_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MEF_MycRas_2DcgPFA_bins15_transition_vectors.csv'))


# Process dfs
dx = [2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5]
binsx = [3, 5, 7, 10, 2, 3, 5, 7, 10, 15, 20, 3, 5, 10, 3, 5]

sig_array = np.zeros((len(dx), 6))
for i in range(len(dx)):
    print('Processing WT...')
    print('D =',dx[i], '|', '# of bins = ', binsx[i])
    print('---')
    wt_sigp, wt_stat = process_samples(wt_df, os.path.join('../data/', celltype, 'nd_pfa'), 'MEF_WT', D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count = 3)
    if np.any(wt_sigp):
        print('Plotting...')
        f = plot_sig_transitions(wt_sigp, os.path.join('../data/', celltype, 'nd_pfa', 'MEF_WT_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '.png'), sigmarkers=True)
    print('Processing MycRas...')
    print('D =',dx[i], '|', '# of bins = ', binsx[i])
    print('---')
    mycras_sigp, mycras_stat = process_samples(mycras_df, os.path.join('../data/', celltype, 'nd_pfa'), 'MEF_MycRas', D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count = 3)
    if np.any(mycras_sigp):
        print('Plotting...')
        f = plot_sig_transitions(mycras_sigp, os.path.join('../data/', celltype, 'nd_pfa', 'MEF_MycRas_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '.png'), sigmarkers=True)

    if np.any(wt_sigp):
        wtSdb = np.any(wt_sigp['reject_H0'])
    else:
        wtSdb = 0

    if np.any(mycras_sigp):
        mycrasSdb = np.any(mycras_sigp['reject_H0'])
    else:
        mycrasSdb = 0

    wtSs = wt_stat['p_val'].min() * 3 # Holm-Bonferroni correction
    mycrasSs = mycras_stat['p_val'].min() * 3 # Holm-Bonferroni correction


    sig_array[i,:] = [wtSdb, mycrasSdb, wtSs, mycrasSs, dx[i], binsx[i]]

sig_df = pd.DataFrame(sig_array)
sig_df.columns = ['wt_sig_db', 'mycras_sig_db', 'wt_p_stationary', 'mycras_p_stationary', 'D', 'n_bins']
sig_df.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'MEF_Aggragate_Sig_ND_PFA.csv'), sep=',', index=False)
