import os
from pfa import *
import pandas as pd
import numpy as np

celltype = 'sims'

D = 2 # number of dimensions to consider
n_bins = 15 # number of bins for course-graining
log_trans = False # set to False for no log transformation
if log_trans:
    log_name = 'logTransformed'
else:
    log_name = ''

# Read in data frames
print('Importing data...')

pwr2rw = read_dataframe(os.path.join('../data/', celltype, 'pwr2urw_split_20.csv'), add_plate=1)
rw2rw = read_dataframe(os.path.join('../data/', celltype, 'rw2rw_split_20.csv'), add_plate=1)

# Plot scatter of 2D transition vectors
pwr2rw_dv = process_dim_vectors(pwr2rw, os.path.join('../data/', celltype, 'nd_pfa', 'pwr2rw_' + str(D) + 'transition_vector_scatter' + log_name + '.png'))
rw2rw_dv = process_dim_vectors(rw2rw, os.path.join('../data/', celltype, 'nd_pfa', 'rw2rw_' + str(D) + 'transition_vector_scatter' + log_name + '.png'))

# Save transition vectors for the 2D, (15,15) bin maps
pwr2rw_sigp, pwr2rw_sv = process_samples(pwr2rw,  os.path.join('../data/', celltype, 'nd_pfa', 'pwr2rw_2DcgPFA_bins15_plotSigT' + '.csv'), D=2, n_bins=15, min_count=3, sv=True)
pwr2rw_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'pwr2rw_2DcgPFA_bins15_transition_vectors.csv'))
rw2rw_sigp, rw2rw_sv = process_samples(rw2rw,  os.path.join('../data/', celltype, 'nd_pfa', 'rw2rw_2DcgPFA_bins15_plotSigT' + '.csv'), D=2, n_bins=15, min_count=3, sv=True)
rw2rw_sv.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'rw2rw_2DcgPFA_bins15_transition_vectors.csv'))

# Process dfs
dx = [2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5]
binsx = [3, 5, 7, 10, 2, 3, 5, 7, 10, 15, 20, 3, 5, 10, 3, 5]

sig_array = np.zeros((len(dx), 4))
for i in range(len(dx)):
    print('---')
    print('Processing FGF2+...')
    print('D =',dx[i], '|', '# of bins = ', binsx[i])
    print('---')
    pwr2rw_sigp = process_samples(pwr2rw, os.path.join('../data/', celltype, 'nd_pfa', 'pwr2rw_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '.csv'), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=3)
    if np.any(pwr2rw_sigp):
        print('Plotting...')
        f = plot_sig_transitions(pwr2rw_sigp,  os.path.join('../data/', celltype, 'nd_pfa', 'pwr2rw_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '.png'), sigmarkers=True)
    print('---')
    print('Processing FGF2-')
    print('D =',dx[i], '|', '# of bins = ', binsx[i])
    print('---')
    rw2rw_sigp = process_samples(rw2rw, os.path.join('../data/', celltype, 'nd_pfa', 'rw2rw_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '.csv'), D = dx[i], n_bins = binsx[i], log_trans=log_trans, min_count=3)
    if np.any(rw2rw_sigp):
        print('Plotting...')
        f = plot_sig_transitions(rw2rw_sigp,  os.path.join('../data/', celltype, 'nd_pfa', 'rw2rw_' + str(dx[i]) + 'DcgPFA_sigTransitions_' + 'bins' + str(binsx[i]).zfill(2) + log_name + '.png'), sigmarkers=True)

    if np.any(pwr2rw_sigp):
        pwr2rwS = np.any(pwr2rw_sigp['reject_H0'])
    else:
        pwr2rwS = 0
    if np.any(rw2rw_sigp):
        rw2rwS = np.any(rw2rw_sigp['reject_H0'])
    else:
        rw2rwS = 0

    sig_array[i,:] = [pwr2rwS,rw2rwS, dx[i], binsx[i]]

sig_df = pd.DataFrame(sig_array)
sig_df.columns = ['pwr2rw_sig', 'rw2rw_sig', 'D', 'n_bins']
sig_df.to_csv(os.path.join('../data/', celltype, 'nd_pfa', 'Sims_Aggragate_Sig_ND_PFA' + log_name + '.csv'), sep=',', index=False)
