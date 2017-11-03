'''
Analyzing real valued transition vectors
'''
import numpy as np
import pandas as pd
import re
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def eval_vector_significance(dim_vectors, D=2, popmean=0):
    '''
    Tests if transition vector magnitudes in each dimension are significantly
    different than a proposed population mean, i.e. 0, using a 1 sample
    t-test

    Parameters
    ----------
    dim_vectors : pandas DataFrame.
        (N * (timesteps - 1)) x (id_vars + 2*D).
        Initial columns are ID variables and a timestep marker, where timesteps
        are integers representing the transition number.
        Next columns are bin locations.
        Final columns are real valued vector components.
    D : integer.
        Number of vector dimensions to test for significance.
    popmean : float.
        Proposed population mean to test against.

    Returns
    -------
    vsig : pandas DataFrame.
        D x 5.
        Contains t statistics and p values for each vector dimension tested.
        ['t_stat', 'p_val', 'vdim', 'v_mean', 'v_sem']
    '''

    vsig = pd.DataFrame(np.zeros((D, 5)))
    vsig.columns = ['t_stat', 'p_val', 'vdim', 'v_mean', 'v_sem']
    for i in range(D):
        v = dim_vectors['vdim' + str(i).zfill(2)]
        t, p = stats.ttest_1samp(v, popmean=popmean)
        vsig.iloc[i,:] = [t, p, i, v.mean(), stats.sem(v)]

    return vsig

def eval_groupwise_vectors(dim_vectors, labels='labels', D = 2):
    '''
    Calculates transition vectors for each unique group in a data set
    i.e. for each cluster, cell type

    Parameters
    ----------
    dim_vectors : pandas DataFrame.
        (N x (id_vars + 2*D).
        Initial columns are ID variables and a timestep marker, where timesteps
        are integers representing the transition number.
        Next columns are real valued red. dim. locations.
        Final columns are real valued vector components.
        Contains a label column with unique labels for each group.
    labels : string.
        name of the label column in dim_vectors.
    D : integer.
        number of vector dimensions to evaluate.
    uniqueT : boolean.
        Calculate vectors for each group at each timepoint.

    Returns
    -------
    group_vectors : pandas DataFrame.
        unique(labels) x (1 + 3*D).
        Contains vector mean, SEM, & pval for each dimension for each group.
        pvals are calculated as the 1 sample ttest against a mean == 0
    '''

    groups = pd.unique(dim_vectors[labels])
    group_vectors = pd.DataFrame( np.zeros( ( len(groups), (1+3*D) ) ) )
    col_names = ['group']
    for d in range(D):
        col_names.append( 'vdim' + str(d).zfill(2) + '_mean')
        col_names.append( 'vdim' + str(d).zfill(2) + '_sem')
        col_names.append( 'vdim' + str(d).zfill(2) + '_p_val')


    group_vectors.columns = col_names

    for i in range(len(groups)):
        g = groups[i]
        g_vectors = dim_vectors[dim_vectors[labels] == g]
        g_vectors.index = range(g_vectors.shape[0])

        g_vsig = eval_vector_significance(g_vectors, D = D, popmean=0.0)

        r = [g]
        for d in range(D):
            r.append(g_vsig.iloc[d, 3])
            r.append(g_vsig.iloc[d, 4])
            r.append(g_vsig.iloc[d, 1])

        group_vectors.iloc[i, :] = r

    return group_vectors
