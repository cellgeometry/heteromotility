#!/usr/bin/python3
'''
Plot Transition vectors for MuSCs atop MuSC/Myoblast PCA
'''
from pfa import *
from trans_vectors import *
import os
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering


# Load MuSC data
# Concat in the same order as the R data
celltype = 'musc'
exp1 = '20160626'
exp2 = '20160701_0'
exp3 = '20160701_1'

fgf = []
nofgf = []
i = 1
print('Importing MuSC data...')
for exp in [exp1, exp2, exp3]:
    fgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'fgf2_exp_motility_statistics_split_20.csv'), add_plate=i))
    nofgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'nofgf2_exp_motility_statistics_split_20.csv'), add_plate=i))
    i += 1

tmp = []
for i in range(len(fgf)):
    tmp.append(fgf[i])
    tmp.append(nofgf[i])

musc_df = pd.concat(tmp)
musc_df.index = range(musc_df.shape[0])

# Load myoblast data
celltype = 'myoblast'
exp1 = '20160623'
exp2 = '20160720'

fgf = []
nofgf = []
print('Importing myoblast data...')
i = 4
for exp in [exp1, exp2]:
    fgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'fgf2_exp_motility_statistics_split_20.csv'), add_plate=i))
    nofgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'nofgf2_exp_motility_statistics_split_20.csv'), add_plate=i))
    i += 1

tmp2 = []
for i in range(len(fgf)):
    tmp2.append(fgf[i])
    tmp2.append(nofgf[i])

myo_df = pd.concat(tmp2)
myo_df.index = range(myo_df.shape[0])

# Combine DFs and obtain transition vectors

comb_df = pd.concat([musc_df, myo_df])
comb_df['cell_type'] = ['MuSC'] * musc_df.shape[0] + ['Myoblast'] * myo_df.shape[0]
comb_df.index = range(comb_df.shape[0])

# Read in full length tracks and define PCA space

celltype = 'musc'
exp1 = '20160626'
exp2 = '20160701_0'
exp3 = '20160701_1'
fgf = []
nofgf = []
i = 1
print('Importing MuSC data...')
for exp in [exp1, exp2, exp3]:
    fgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'fgf2_exp_motility_statistics.csv'), add_plate=i))
    nofgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'nofgf2_exp_motility_statistics.csv'), add_plate=i))
    i += 1

flmusc_df = pd.concat( [pd.concat(fgf), pd.concat(nofgf)] )

celltype = 'myoblast'
exp1 = '20160623'
exp2 = '20160720'

fgf = []
nofgf = []
print('Importing myoblast data...')
i = 4
for exp in [exp1, exp2]:
    fgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'fgf2_exp_motility_statistics.csv'), add_plate=i))
    nofgf.append(read_dataframe(os.path.join('../data/', celltype, exp, 'nofgf2_exp_motility_statistics.csv'), add_plate=i))
    i += 1

flmyo_df = pd.concat( [pd.concat(fgf), pd.concat(nofgf)] )

flength_df = pd.concat([flmusc_df, flmyo_df])
flength_df_dr, pca = reduce_dims(flength_df.iloc[:,3:57], D=2, eigenv=True)

# Read in labels exported from R clustering
labels = pd.read_csv('../data/comparisons/muscmyo/muscmyo_hclust_labels.csv')
# apply labels to comb_df for export
comb_df_id = split_col(comb_df.iloc[:,:3], col_name='cell_id')

comb_df_dr = pd.DataFrame( reduce_dims(comb_df.iloc[:,3:57], w=pca) )
comb_df_dr.columns = ['PC1', 'PC2']

comb_df_drid = pd.concat([comb_df_id, comb_df_dr], 1)

tmp3 = []
for i in range(np.array(comb_df_drid.timestep, dtype='int32').max()+1):
    tv = comb_df_drid[comb_df_drid.timestep == str(i)]
    tv.index = range(tv.shape[0])
    tv.insert(0, 'labels', labels)
    tmp3.append(tv)
lab_comb_df = pd.concat(tmp3)
lab_comb_df.to_csv('../data/comparisons/muscmyo/muscmyo_labeled_state_df.csv', index = False)


# Calculate transition vectors in the PC space defined using full length tracks

print('Calculating transition vectors...')
dim_vectors = process_dim_vectors(comb_df, os.path.join('../data/comparisons/muscmyo', 'MuSC_Myo_Vector_Scatter.png'), n_ids=3, w=pca)

tmp4 = []
for i in range(np.array(dim_vectors.timestep, dtype='int32').max()+1):
    tv = dim_vectors[dim_vectors.timestep == str(i)]
    tv.index = range(tv.shape[0])
    tv.insert(0, 'labels', labels)
    tmp4.append(tv)

dim_vectors = pd.concat(tmp4)

# Separate MuSC and Myoblast red. dim. vectors
musc_dv = dim_vectors.iloc[:int(musc_df.shape[0] * (2/3)),:]
myo_dv = dim_vectors.iloc[int(musc_df.shape[0] * (2/3)):, :]

# Test for transition significance for population
musc_vsig = eval_vector_significance(musc_dv, D = 2)
myo_vsig = eval_vector_significance(myo_dv, D = 2)

# Check transition significance for subpopulations
group_vectors_all = eval_groupwise_vectors(dim_vectors, D=2)
group_vectors_t0 = eval_groupwise_vectors(dim_vectors[dim_vectors.timestep=='0'], D=2)
group_vectors_t1 = eval_groupwise_vectors(dim_vectors[dim_vectors.timestep=='1'], D=2)

group_vectors_all.insert(0, 'timestep', 'all')
group_vectors_t0.insert(0, 'timestep', 0)
group_vectors_t1.insert(0, 'timestep', 1)

group_vectors = pd.concat([group_vectors_all, group_vectors_t0, group_vectors_t1])
group_vectors.to_csv(os.path.join('../data/comparisons/muscmyo/', 'muscmyo_groupwise_transition_vectors.csv'), index = False)
