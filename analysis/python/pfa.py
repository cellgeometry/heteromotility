'''
Probability Flux Analysis

Implements PFA in N-dimensions for a given feature matrix

Workflow:

DataFrame >> DimReduction >> N-dim binning

for cell in DataFrame:
    For window in timecourse:
        identify initial bin
        determine next bin
        record transition vector

Output:

Cell_i, Vector_0, Vector_1, .. Vector_N
'''

import numpy as np
import pandas as pd
import re
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import math
from statsmodels.stats import multitest
import itertools
import os
import sys
from kaufmann import Kaufmann2003Solve

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')
sns.set_palette('muted')

def read_dataframe(df_location, scale=True, add_plate=False):
    '''
    Reads feature data CSV into DataFrame

    Parameters
    ----------
    df_location: string
        path to feature data CSV.
    scale : boolean.
        perform centering and scaling by removing feature means and scaling
        to unit variance.

    Returns
    -------
    df : pandas DataFrame.
        N x M feature data frame.
    '''

    df = pd.read_csv(df_location)
    if add_plate:
        df.insert(0, 'plate', add_plate)

    return df

def reduce_dims(df, D = 3, scale = True, log_trans = False, eigenv=False, w=None):
    '''
    Reduces dimensionality of an N x M feature DataFrame

    Parameters
    ----------
    df : pandas DataFrame.
        N x M feature data frame.
    D : integer.
        number of dimensions from reduced dimensional space to keep.
    log_trans : boolean.
        Perform a log transform on reduced dimensional features.
    eigenv : boolean.
        Return PCA object with eigenvector matrix.
    w : PCA object.
        sklearn.decomposition.PCA object with fitted eigenvector matrix.
        Used to transform new points into a previously defined space.

    Returns
    -------
    df_dr : ndarray.
        N x D reduced dimensional feature data frame, where D is a parameter.
    '''
    if scale and not log_trans:
        scaler = StandardScaler()
        df = scaler.fit_transform(df)

    if w:
        df_dr = w.transform(df)
    else:
        pca = PCA(n_components = D)
        pca.fit(df)
        df_dr = pca.fit_transform(df)
    if log_trans:
        df_dr += np.abs(df_dr.min()) + 0.001
        df_dr = np.log(df_dr)

    if eigenv:
        return df_dr, pca
    else:
        return df_dr

def split_col(df, col_name, delimiter='-'):
    '''
    Splits a specified column of a dataframe around
    a provided delimiter.
    '''
    s = pd.Series(df[col_name])
    split = s.str.split(pat=delimiter)
    cell = split.str.get(0)
    timestep = split.str.get(1)

    df_split = df.drop(col_name, 1)
    df_split['cell'] = cell
    df_split['timestep'] = timestep
    return df_split

def bin_dimension(values, n_bins, sigfig=4):
    '''
    Determines which course-grained bin along a dimension
    values occupy.

    Parameters
    ----------
    values : ndarray. N x 1.
        Values to determine bin coordinates.
    n_bins : integer.
        number of bins for course-graining the dimension.
    sigfig : integer.
        number of significant figures to consider.

    Returns
    -------
    bins : ndarray. N x 1.
        Corresponding bin locations for each value.

    Notes
    -----
    * Rounds all values to four sigfigs.
    * If a value is on the border of two bins, randomly chooses one location.
    '''

    values = np.round(np.reshape(values, (values.shape[0], 1)), sigfig) # round to 4 sigfigs
    dim_range = np.array([ np.min(values), np.max(values) ])
    bin_sz = np.abs((dim_range[1] - dim_range[0])) / n_bins
    # Generate a matrix N x n_bins
    # Values of each column are the min value of a given bin
    bin_mat = np.zeros([values.shape[0], n_bins])
    for i in range(n_bins):
        bin_mat[:,i] = dim_range[0] + [i * bin_sz]
    # Generate a matrix of values, N x n_bins
    val_mat = np.tile(values, (1, n_bins))
    # Compare values elementwise to determine which bin a value falls into
    gt_mat = (val_mat >= np.round(bin_mat, sigfig))
    # Compare to max bin values to generate the corresponding binary matrix
    lt_mat = ( val_mat <= np.round((bin_mat+bin_sz), sigfig) )
    # Overlap between these matrices is the bin location of a value
    # store values in a course-grained matrix
    cg_mat = (gt_mat == lt_mat)
    # Check for values on the border of a bin, which will return
    # 2 True indices in cg_mat
    double_count = np.where( np.sum(cg_mat, axis=1) == 2 )
    double_count = double_count[0] # pop array from two val tuple
    # For each double_counted value, randomly choose a bin to select
    select = np.round(np.random.random()).astype('int32')
    for i in range(double_count.shape[0]):
        true_idx = np.where(cg_mat[double_count[i], :])[0]
        # Set one of the two vals to False
        cg_mat[double_count[i], np.int32(true_idx[select])] = False

    idx, bin_loc = np.where(cg_mat)
    return idx, bin_loc

def gen_state_locations(df_id, df_dr, n_bins=15):
    '''
    Records the location of each sample at each timestep
    by binning low dimensional space into course-grained states

    Parameters
    ----------
    df_id : pandas DataFrame.
        Contains unique ID variables, with the final column indicating the time
        window of the observation, recorded as integers.
    df_dr : pandas DataFrame or ndarray.
        Reduced dimensional N x D feature matrix.
    n_bins : integer.
        Number of equally sized bins to course grain each dimension.

    Returns
    -------
    state_locations : pandas DataFrame.
        N x (id_vars + D) DataFrame, where the initial columns are all unique
        ID variables and the trailing columns are locations in binned reduced
        dimensional space.
        Final ID column is the timestep, named 'timestep'.
    dim_locations : pandas DataFrame.
        N x (id_vars + D) DataFrame, where the initial columns are all unique
        ID variables and the trailing columns are real valued locations in
        reduced dimensional space.
    '''

    time_name = df_id.columns[-1]
    tmp_df = pd.DataFrame(np.zeros( [df_dr.shape[0], df_dr.shape[1]] ))
    dim_names = []
    for i in range(df_dr.shape[1]):
        dim_names.append('dim' + str(i).zfill(2))
    tmp_df.columns = dim_names
    state_locations = pd.concat([df_id, tmp_df], axis=1)

    for d in range(df_dr.shape[1]):
        idx, bin_loc = bin_dimension(df_dr[:,d], n_bins=n_bins)
        state_locations.ix[:,(df_id.shape[1] + d)] = bin_loc

    tmp2_df = pd.DataFrame(df_dr)
    tmp2_df.columns = dim_names
    dim_locations = pd.concat([df_id, tmp2_df], 1)

    return state_locations, dim_locations

# Check if the system is stationary
# i.e. if proportion of samples in each state is consistent over timesteps

def nCr(n, r):
    '''
    Combinations of size r in sample of size n
    '''
    f = math.factorial
    return f(n) / ( f(r)*f(n-r) )

def check_stationary(state_locations, D=3):
    '''
    Determines if a system is stationary over time by using a chi2 test
    for difference between the contigency tables of state distribution at
    each time point in the series.

    Parameters
    ----------
    state_locations : pandas DataFrame.
        N x (id_vars + D) DataFrame, where the initial columns are all unique
        ID variables and the trailing columns are locations in binned reduced
        dimensional space.
        Final ID column is the timestep, named 'timestep'.
    D : integer.
        Dimensionality of state space.
    occupied_only : boolean.
        removes bins that have no occupants at either timestep from the
        contingency table before performing a chi2 test.

    Returns
    -------
    stationary_stats : pandas DataFrame.
        C(T, 2) x 4 DataFrame, where T is the maximum timestep.
        Each row contains statistics for differences between each pair of intervals.
        Columns 1, 2 denote the initial and subsequent timestep being compared.
        Column 3, 4 denote the Chi2 statistic and p-value.
    '''

    T = np.int32( state_locations['timestep'].max() )
    comparisons = list(itertools.combinations(range(T+1), 2))

    stationary_stats = pd.DataFrame(np.zeros((len(comparisons), 4)))
    stationary_stats.columns = ['t0', 'tT', 'chi2', 'p_val']


    k = 0
    for t0, t1 in comparisons:
        # Build contigency table
        # n_bins**D x D+2 DataFrame, where first D cols are state bins
        # and last two cols are t0 cell count and tT cell count
        n_bins = np.int32(state_locations['dim00'].max())

        cont_table = pd.DataFrame(np.zeros((n_bins**D, D+2)))
        col_names = list(state_locations.columns[-D:])
        col_names += ['t0', 'tT']
        cont_table.columns = col_names

        # Find all possible combinations of n_bins with itertools
        A = []
        for d in range(D):
            A.append(list(range(n_bins)))

        p = itertools.product(*A) # uses the splat operator to expand the list
        cont_table.iloc[:,:D] = np.array(list(p))

        for r in range(cont_table.shape[0]):
            sl_t0 = state_locations.iloc[np.int32(state_locations['timestep']) == t0, :]
            sl_t1 = state_locations.iloc[np.int32(state_locations['timestep']) == t1, :]
            count_t0 = np.sum( np.sum(sl_t0.iloc[:,-D:] == cont_table.iloc[r,:D], 1) == D )
            count_t1 = np.sum( np.sum(sl_t1.iloc[:,-D:] == cont_table.iloc[r,:D], 1) == D )

            cont_table.iloc[r,-2:] = np.array([count_t0, count_t1])

        # Test contingency table for independence of timepoints
        rc = np.array(cont_table.iloc[:,-2:])
        rc = rc[np.sum(rc, 1)>0,:]
        chi2, p, dof, expected = stats.chi2_contingency(rc.T)

        stationary_stats.iloc[k,:] = [t0, t1, chi2, p]
        k += 1

    return stationary_stats

def gen_state_vectors(state_locations, D=3):
    '''
    Calculates transition vectors between timesteps for each sample.

    Parameters
    ----------
    state_locations : pandas DataFrame.
        (N*timesteps) x (id_vars + D) DataFrame, where the initial columns are all unique
        ID variables and the trailing columns are locations in binned reduced
        dimensional space.
        The final ID variable column is assumed to be the timestep, named 'timestep'.
        ie. ID1, ID2,.., IDn, Time, Dim1, Dim2, ..., DimN
    D : integer.
        Number of feature dimensions for consideration.


    Returns
    -------
    state_vectors : pandas DataFrame.
        (N * (timesteps - 1)) x (id_vats + 2*D).
        Initial columns are ID variables and a timestep marker, where timesteps
        are integers representing the transition number.
        Next columns are bin locations.
        Final columns are vector components.
    '''
    # Sort values in state locations based on id_vars
    nb_ids = len(state_locations.columns) - D
    id_names = state_locations.columns[:nb_ids]
    vdim_names = []
    for i in range(D):
        vdim_names.append('vdim' + str(i).zfill(2))
    state_locations = state_locations.sort_values(list(id_names))
    state_locations.index = range(state_locations.shape[0])
    # Max timestep = (how many timesteps per cell - 1)
    T = np.max(np.int32(state_locations['timestep']))
    N = int(state_locations.shape[0]/(T+1))

    state_vectors = pd.DataFrame()

    for i in range(0, state_locations.shape[0], T+1):
        cell = np.array(state_locations.iloc[i:(i+T+1),nb_ids:])
        # Shift array by one timestep to find vectors
        vec = cell[1:(T+1),:] - cell[0:T,:]
        vec_df = pd.DataFrame(vec)
        vec_df.columns = vdim_names
        id_df = state_locations.iloc[i:(i+T),:]
        id_df.index = range(T)
        tmp_df = pd.concat([id_df, vec_df], 1)
        state_vectors = state_vectors.append(tmp_df)

    state_vectors.index = range(state_vectors.shape[0])
    return state_vectors

# Check pairwise state detailed balance breaking

def pairwise_balance(state_vectors, n_bins = 15, D = 3, n_ids=3, min_count = 6):
    '''
    Determine if pairwise transition rates obey detailed balance
    using the binomial test to check 1:1 ratio of pairwise transitions

    Parameters
    ----------
    state_vectors : pandas DataFrame.
        (N * (timesteps - 1)) x (id_vats + 2*D).
        Initial columns are ID variables and a timestep marker, where timesteps
        are integers representing the transition number.
        Next columns are bin locations.
        Final columns are vector components.
    n_bins : integer.
        number of bins used for each dimension.
    D : integer.
        number of dimensions considered.
    n_ids : integer.
        number of ID variables, excluding the timestep.
    min_count : integer.
        minimum number of transitions in a pair to test by binomial test.
        Defualt = 6. 6 trials is minimum for sufficient power at alpha = 0.05.

    Returns
    -------
    p_vals : pandas DataFrame.
        P x (1+ 2*D) array of P values for a binomial test of balance between
        transitions.
        Col 0 is the p value. Subsequent columns are the initial position and
        position after transition.

    Notes
    -----
    Constructs a state transition space hypershape of dimensionality D*2
    The first D dimensions represent a samples initial state location,
    the second D dimensions represent the location after it's transition
    The value in each unit of the space records how many samples underwent a
    transition between those locations.
    '''

    # Create a hyperspace to record all possible transitions
    transition_space = np.zeros(np.tile(n_bins, D*2))
    trans_points = np.max(np.int32(state_vectors['timestep'])) + 1

    for i in range(0, state_vectors.shape[0], trans_points):
        cell = state_vectors.iloc[i:(i+trans_points),(n_ids+1):]
        for j in range(trans_points):
            p0 = np.array(cell.iloc[j][:D]).astype('int32')
            v0 = np.array(cell.iloc[j][D:]).astype('int32')
            p1 = p0 + v0
            idx = tuple(np.concatenate([p0, p1]))
            transition_space[idx] += 1

    # Check to ensure some transitions are above the minimum criteria
    if transition_space.max() <= min_count:
        print('No transitions met minimum significance criteria.')
        print('Consider using fewer bins or dimensions.')
        return pd.DataFrame(np.zeros((0,0)))

    # specific transitions are rows, columns are each dims index
    # only consider transitions with a meaningful number of samples
    # i.e. with enough power for a binomial test
    trans_F = np.array(np.where(transition_space > min_count)).T
    # find the reverse transitions
    trans_R = np.hstack([ trans_F[:,D:], trans_F[:,:D] ])
    # Find unique transitions so we don't test the same ones twice
    trans_both = np.vstack([trans_F, trans_R])
    trans_U = np.vstack(set(map(tuple, trans_both)))

    p_vals = np.zeros([trans_U.shape[0], (1 + D*2 + 2)])

    for i in range(trans_U.shape[0]):
        idx = tuple(trans_U[i,:])
        rev_idx = idx[D:] + idx[:D]
        t0 = transition_space[idx]
        t1 = transition_space[rev_idx]
        p = stats.binom_test(t0, np.sum([t0, t1]))
        p_vals[i,:] = tuple([p]) + idx + (t0, t1)

    p_vals = pd.DataFrame(p_vals)
    col_names = ['pval']
    for i in range(D):
        col_names.append( 't0_dim'+str(i).zfill(i) )
    for i in range(D):
        col_names.append( 't1_dim'+str(i).zfill(i) )
    col_names.append('state0')
    col_names.append('state1')
    p_vals.columns = col_names
    return p_vals.sort_values('pval')

def fdr_correction(p_vals, alpha = 0.05):
    '''
    Corrects p-values using the Benjamini-Hochberg method to a desired
    false discovery rate (FDR).

    Parameters
    ----------
    p_vals : pandas DataFrame.
        P x (1+D) array of P values for a binomial test of balance between
        transitions.
        Col 0 is the p value. Subsequent columns are the initial position and
        position after transition.
    alpha : float, (0, 1).
        alpha assurance value to control FDR.

    Returns
    -------
    sig_p_vals : pandas DataFrame.
        sigP x (2 + 2*D + 2) array of corrected p-values, controlling the FDR at the
        specified alpha.
        First two columns the rejection of null hypothesis (boolean), and the
        corrected pval.
        Following columns are the location of initial state and destination
        state, each with D coordinates.
        Final two columns are the number of cells making each transition.comb_df_dr, pca = reduce_dims(comb_df.iloc[:,3:57], D = 2, eigenv=True)

    '''

    p = np.asarray(p_vals['pval'])
    rejected, pval_corrected = multitest.fdrcorrection(p, alpha = alpha, method = 'indep', is_sorted = True)
    tmp_df = pd.DataFrame([rejected, pval_corrected]).T
    tmp_df.columns = ['reject_H0', 'pval_corr']
    p_vals.index = range(p_vals.shape[0])
    sig_p_vals = pd.concat([tmp_df, p_vals.iloc[:,1:]], 1) # join by cols

    return sig_p_vals


def process_samples(df, output_location, exp_name,
                    n_ids=3, D=3, col_name='cell_id',
                    n_bins=15, min_count=6, max_feat=57,
                    alpha=0.05, log_trans=False,
                    verbose=True, sv=False):
    '''
    Process a data frame to find significantly unbalanced state transitions
    by N-dimensional course-grained probability flux analysis (Nd-cgPFA).

    Parameters
    ----------
    df : pandas DataFrame.
        N x M dataframe, where the first n_ids columns are unique identifiers,
        with the final identifier col_name containing sample ID and timestcomb_df_dr, pca = reduce_dims(comb_df.iloc[:,3:57], D = 2, eigenv=True
        information, separated by a '-'.
        Remaining M-n_ids columns are unique features.
    output_location : string.
        path for CSV output.
    exp_name : string.
        Prefix for CSV file output.
    n_ids : integer.
        number of unique identifiers in df.
    D : integer.
        number of dimensions to consider for N-dim PFA.
    col_name : string.
        name of column containing sample ID and timestep info, separated by '-'.
    n_bins : integer.
        number of bins for course-graining each dimension.
    min_count : integer.
        minimum number of samples transitioning between two bins to consider
        for significance testing.
        Setting to 1 tests all bins with any transitions.
        Values >1 ignore transitions with too few cells to reach signficance.
    max_feat : integer.
        maximum feature number in the raw df for consideration.
        Setting to None uses all features.
    alpha : float, (0, 1).
        Assurance level to control FDR during multiple hypothesis correction.
    log_trans : boolean.
        Perform a log transform on reduced dimensional features before finding
        transitions.
    sv : boolean.
        return state_vectors.
    verbose : boolean.
        Print signficance or non-significance to stdout.

    Returns
    -------
    sig_p_vals
    state_vectors
    '''
    df_dr = reduce_dims(df.iloc[:,n_ids:max_feat], D = D, log_trans = log_trans)
    df_id = split_col(df.iloc[:,:n_ids], col_name=col_name)
    state_locations, dim_locations = gen_state_locations(df_id, df_dr, n_bins=n_bins)

    # Check stationary nature of system
    stationary_stats = check_stationary(state_locations, D = D)
    stationary_stats.to_csv(os.path.join(output_location, exp_name + str(D)+'_D_' + 'bins' + str(n_bins).zfill(2) + '_stationary_stats.csv'), sep=',', index = False)

    state_vectors = gen_state_vectors(state_locations, D = D)
    p_vals = pairwise_balance(state_vectors, n_bins = n_bins, D = D, min_count = min_count)
    # check if no transitions met min_count criteria
    if not np.any(p_vals) and sv==False:
        return None, stationary_stats
    elif not np.any(p_vals) and sv==True:
        return None, state_vectors
    sig_p_vals = fdr_correction(p_vals, alpha = alpha)
    if verbose and np.any(sig_p_vals['reject_H0']):
        print('Significant detailed balance breaking.')
    elif verbose and not np.any(sig_p_vals['reject_H0']):
        print('No significant detailed balance breaking.')
    sig_p_vals.to_csv(os.path.join(output_location, exp_name + str(D) + 'DcgPFA_sigTransitions_' + 'bins' + str(n_bins).zfill(2) + '.csv'), sep=',', index = False)



    if sv:
        return sig_p_vals, state_vectors
    else:
        return sig_p_vals, stationary_stats

def process_stationary(df, output_location, exp_name,
                    n_ids=3, D=3, col_name='cell_id',
                    n_bins=15, max_feat=57):
    '''
    Process a data frame to check if transitions are stationary.

    Parameters
    ----------
    df : pandas DataFrame.
        N x M dataframe, where the first n_ids columns are unique identifiers,
        with the final identifier col_name containing sample ID and timestcomb_df_dr, pca = reduce_dims(comb_df.iloc[:,3:57], D = 2, eigenv=True
        information, separated by a '-'.
        Remaining M-n_ids columns are unique features.
    output_location : string.
        path for CSV output.
    exp_name : string.
        Prefix for CSV file output.
    n_ids : integer.
        number of unique identifiers in df.
    D : integer.
        number of dimensions to consider for N-dim PFA.
    col_name : string.
        name of column containing sample ID and timestep info, separated by '-'.
    n_bins : integer.
        number of bins for course-graining each dimension.

    max_feat : integer.
        maximum feature number in the raw df for consideration.
        Setting to None uses all features.


    Returns
    -------
    sig_p_vals
    state_vectors
    '''
    df_dr = reduce_dims(df.iloc[:,n_ids:max_feat], D = D)
    df_id = split_col(df.iloc[:,:n_ids], col_name=col_name)
    state_locations, dim_locations = gen_state_locations(df_id, df_dr, n_bins=n_bins)

    # Check stationary nature of system
    stationary_stats = check_stationary(state_locations, D = D, occupied_only=True)
    stationary_stats.to_csv(os.path.join(output_location, exp_name + str(D)+'_D_' + 'bins' + str(n_bins).zfill(2) + '_stationary_stats.csv'), sep=',', index = False)

    return stationary_stats


#----------------------
# Dwell Time Functions
#----------------------

def check_dwell_time(cell):
    '''
    Determines the dwell time of a given cell in a state. Checks to see if
    occupancies are concurrent, or if they are split.
    Returns a list of dwell times.

    Parameters
    ----------
    cell : pandas DataFrame.
        A single cell in the format of the `state_locations` DataFrame.
        i.e., the following with only a single value in cell

        (N*timesteps) x (id_vars + D) DataFrame, where the initial columns are all unique
        ID variables and the trailing columns are locations in binned reduced
        dimensional space.
        The final ID variable column is assumed to be the timestep, named 'timestep'.
        ie. ID1, ID2,.., IDn, Time, Dim1, Dim2, ..., DimN

    Returns
    -------
    dwell : ndarray.
        array of integers of discrete time units the cell dwelled in the state.
    '''

    minT = int(cell['timestep'].min())
    maxT = int(cell['timestep'].max())
    # check if steps occured concurrently
    steps = np.array(cell['timestep'].sort_values().astype('int'))
    concurrent = np.all( steps == range(minT, maxT+1) )

    if concurrent:
        dwell = [int(cell.shape[0])]
    else:
        # check where the range has breaks
        # todo : This is O(n^2) and could be improved
        breaks = [x for x in range(minT, maxT+1) if x not in steps]

        dwell = []
        for i in range(len(breaks)):
            if i == 0:
                dwell.append( np.sum(steps < breaks[i]) )
            elif i == len(breaks)-1:
                dwell.append( np.sum(steps > breaks[i]) )
            else:
                b0 = breaks[i]
                b1 = breaks[i+1]
                dwell.append( np.sum( np.logical_and(steps > b0, steps < b1) ) )
                dwell = [x for x in dwell if x > 0] # remove 0's

    return np.array(dwell)

def gen_dwell_times(state_locations, D = 3):
    '''
    Generates a distribution of dwell times for each state

    Parameters
    ----------
    state_locations : pandas DataFrame.
        (N*timesteps) x (id_vars + D) DataFrame, where the initial columns are all unique
        ID variables and the trailing columns are locations in binned reduced
        dimensional space.
        The final ID variable column is assumed to be the timestep, named 'timestep'.
        ie. ID1, ID2,.., IDn, Time, Dim1, Dim2, ..., DimN

    Returns
    -------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    '''

    # Identify occupied states in the space
    id_vars = state_locations.shape[1] - D
    all_states = np.array( state_locations.iloc[:, -D:] )
    occupied = np.unique( all_states, axis=0 )

    # Generate a dwell times DataFrame with D columns identifying occupied bins
    # and T columns identifying the discrete dwell time of occupants

    T = int(state_locations['timestep'].max()) + 1
    time_bins = np.zeros([occupied.shape[0], T])
    dwell_times = pd.DataFrame(np.hstack([occupied, time_bins]))
    time_bins_names = []
    for i in range(T):
        time_bins_names.append('t' + str(i).zfill(2))
    dwell_times.columns = list(state_locations.columns[-D:]) + time_bins_names

    # For each state, consider each cell that occupies it and determine how long
    # the cell dwelled in that state

    for i in range(occupied.shape[0]):
        state = occupied[i,:] # state as a vector
        bool_idx = np.all(np.equal(all_states, state), axis=1)
        occupants = state_locations.iloc[bool_idx,:]
        # sort the occupants to get unique cells in row order
        occupants = occupants.sort_values(['plate', 'Well/XY', 'cell', 'timestep'])
        # add D+1 to drop timestep col
        unique_occupants = occupants.iloc[:,:-(D+1)].drop_duplicates(occupants.iloc[:,:-(D+1)])

        for j in range(unique_occupants.shape[0]):
            bool_idx = (occupants.iloc[:,:-(D+1)] == unique_occupants.iloc[j,:]).all(1)
            cell = occupants[bool_idx]
            dwell = check_dwell_time(cell)
            for d in dwell:
                dwell_times.iloc[i, (d-1 + D)] = dwell_times.iloc[i, (d-1 + D)] + 1

    return dwell_times

def calc_mean_dwell_time(dwell_times, T):
    '''
    Calculates the mean dwell times for each state.

    Parameters
    ----------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    T : integer.
        number of timesteps.

    Returns
    -------
    mean_dwell_times : pandas Series.
        mean dwell time for each state, indexed per `dwell_times`.
    '''

    times = dwell_times.iloc[:,-T:]
    mean_dwell_times = times.multiply(np.arange(1,T+1), axis=1).sum(1) / times.sum(1)

    return mean_dwell_times

def calc_dwell_time_constants(dwell_times, T, N=1, overall_tc=False):
    '''
    Calculates time constants for state decay by fitting exponential sums
    'observed occupants', as determined by the cumulative sum of dwell time
    observations.

    i.e.) If you look at state S_i at t=1, you will see cells with dwell
    d = {1, 2, 3}.

    Uses the numerical method described in :
        Kaufmann 2003, arXiv:physics/0305019
        Fitting a sum of exponentials to numerical data

    With implementation by:
        Kyle Hurston (@khurston)
        https://github.com/khuston/pyfitdecay

    Fits a sum of exponentials of the form:

        y(t) = a_0 + sum_i^N( a_i * exp((-b_i) t) )

    such that

        tau_i = b_i**-1

    in the standard form for exponential decay functions

        y(t) = C + sum_i^N( A*exp(-t/tau_i) )

    Parameters
    ----------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    T : integer.
        number of timesteps.
    N : integer.
        number of exponentials to fit per state.
    overall_tc : boolean.
        return a value for the time constant of decay across all states.

    Returns
    -------
    time_constants : ndarray.
        NumberOfStates x N array of time constants tau.
    tau_overall : float, optional.
        time constant for decay across all states.
        only returned if overall_tc == True.
    '''

    time_constants = np.zeros((dwell_times.shape[0], N))

    x = np.arange(1, T+1)
    for i in range(dwell_times.shape[0]):
        d = np.array(dwell_times.iloc[i,-T:])
        # Set y as the # of cells observed in the state for a given dwell time
        # by taking the reverse cumsum -- i.e. count dwells of 2 in observation
        # at time 1
        y = np.cumsum(d[::-1])[::-1]
        # Fit a sum of exponentials with Kauffman method
        a, b = Kaufmann2003Solve(N=N, x=x, y=y)
        time_constants[i, :] = b

    if overall_tc:
        a, b = Kaufmann2003Solve(N=N, x=x, y=dwell_times.iloc[:,-T:].sum(0))
        return time_constants**-1, b**-1
    else:
        return time_constants**-1

def expd(x, a0, a1, b):
    '''
    Return value for exponential distribution with parameters
    `a0`, `a1`, `b`, at `x`
    '''
    return a0 + a1*np.exp(-b*x)

def evol_dwell_time_constants(dwell_times, T, overall_tc=True):
    '''
    Fits a single exponential of the form

            y(t) = sum_i^N( a_i * exp((-tau_i) t) )

    using curve fitting.

    Parameters
    ----------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    T : integer.
        number of timesteps.
    overall_tc : boolean.
        return a value for the time constant of decay across all states.

    Returns
    -------
    time_constants : ndarray.
        NumberOfStates x 1 array of time constants tau.
    scales : ndarray.
        scale parameters, NumberOfStates x 1.
    tau_overall : float, optional.
        time constant for decay across all states.
        only returned if overall_tc == True.
    scale_overall : float, optional.
    '''
    from diffevol import ExpFitDiffEvol
    time_constants = np.zeros(dwell_times.shape[0])
    scales = np.zeros(dwell_times.shape[0])
    x = np.arange(1, T+1)
    for i in range(dwell_times.shape[0]):
        d = np.array(dwell_times.iloc[i,-T:])
        # Set y as the # of cells observed in the state for a given dwell time
        # by taking the reverse cumsum -- i.e. count dwells of 2 in observation
        # at time 1
        y = np.cumsum(d[::-1])[::-1]
        a, tau = ExpFitDiffEvol(1, x, y)

        time_constants[i] = tau
        scales[i] = a
    if overall_tc:
        a, tau = ExpFitDiffEvol(1, x, dwell_times.iloc[:,-T:].sum(0))
        return time_constants**-1, scales, tau**-1, a
    else:
        return time_constants**-1, scales

def test_dwell_lilliefors(dwell_times, T = 3):
    '''
    Tests if a distribution is exponentially distributed using Lilliefors' test.
    p-values < alpha_level indicate the distribution is exponential.

    Parameters
    ----------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    T : integer.
        number of timesteps.

    Returns
    -------
    d_ks : float.
        chi^2 goodness of fit statistic.
    pvalue : float.
        p-value returned by K-S test.
    '''

    from lilliefors import lilliefors

    dwell_sums = np.array(dwell_times.iloc[:,-T:].sum(0))
    dwells = []
    for i in range(dwell_sums.shape[0]):
        dwells.append(np.repeat(i+1, dwell_sums[i]))

    x = np.concatenate(dwells)
    d_ks, p = lilliefors(x, dist='exp')

    return d_ks, p

def test_dwell_chisq(dwell_times, T = 3):
    '''
    Test dwell time distribution against exponential theoretical distribution
    using Pearson's chi^2 test of the *source x bin* contingency table.

    Parameters
    ----------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    T : integer.
        number of timesteps.

    Returns
    -------
    chi2 : float.
        chi2 test statistic.
    pval : float.
        p-value of chi2 test.
    theoretical : ndarray.
        theoretical distribution of dwell times.
    '''

    dwell_sums = np.array(dwell_times.iloc[:,-T:].sum(0))
    dwells = []
    for i in range(dwell_sums.shape[0]):
        dwells.append(np.repeat(i+1, dwell_sums[i]))
    # Generate binned theoretical values
    x = np.concatenate(dwells)
    ex = stats.expon(scale=1./np.mean(x))
    exp_var = ex.rvs(size=np.int(1e6))
    theoretical = []
    for i in range(T):
        theoretical.append( np.logical_and(exp_var >= i, exp_var < i + 1).sum() )
    theoretical = (np.array(theoretical) * x.size/1e6).astype('int32')
    # final bin should include all vals > T, as these cells would be observed to
    # dwell for 'only' T units
    theoretical[T-1] += x.size - theoretical.sum()

    observed = np.vstack([dwell_sums, theoretical])
    chi2, pval, _, _ = stats.chi2_contingency(observed)

    return chi2, pval, theoretical

def scale_based_theor_exp(dwell_times, tc, scale, T = 3):
    '''
    Test dwell time distribution against exponential theoretical distribution
    using Pearson's chi^2 test of the *source x bin* contingency table.

    Parameters
    ----------
    dwell_times : pandas DataFrame.
        N x (D + T), where N is the number of bins and D is the number of
        dimensions, and T is the number of discrete time steps
        First D columns identify the state, following T columns represent the
        number of cells with the corresponding dwell time in that state
        Effectively, encodes the bins of a discrete time histogram.
    scale : float.
        scale parameter for expon theoretical curve.
    T : integer.
        number of timesteps.

    Returns
    -------
    chi2 : float.
        chi2 test statistic.
    pval : float.
        p-value of chi2 test.
    theoretical : ndarray.
        theoretical distribution of dwell times.
    '''
    dwell_sums = np.array(dwell_times.iloc[:,-T:].sum(0))
    dwells = []
    for i in range(dwell_sums.shape[0]):
        dwells.append(np.repeat(i+1, dwell_sums[i]))

    def theor_ed(x, scale, tc):
        tau = tc**-1
        return scale*np.exp(-x*tau)

    for x in range(1, T+1):
        theoretical += [theor_ed(x, scale, tc)]
    theoretical = np.array(theoretical).astype('int32')



def plot_dwell_v_occupancy(mean_dwell, sum_dwell, output_location, width=4, height=3):
    g = sns.jointplot(mean_dwell, sum_dwell)
    g.set_axis_labels(xlabel="Mean Dwell Time", ylabel="Total Dwelling Cells")
    plt.savefig(output_location)
    plt.close()
    return

def plot_tc_v_occupancy(time_constants, sum_dwell, output_location, width=4, height=3):
    '''
    Note: can only been used when N=1 exponentials are fit to dwell time data
    '''
    g = sns.jointplot(time_constants.squeeze(), sum_dwell)
    g.set_axis_labels(xlabel="Time Constant", ylabel="Total Dwelling Cells")
    plt.savefig(output_location)
    plt.close()
    return

def plot_dwell_times(dwells, theor, times, output_location, width=4, height=3):
    '''
    Parameters
    ----------
    dwells : array-like of cells per time unit.
    theor : array-like of theoretical dwells per time unit.
    times : array-like of time units.
    '''

    timesr = np.tile(times, 2)
    counts = np.concatenate([dwells, theor])
    source = ['Observed']*dwells.size + ['Theoretical']*theor.size
    df = pd.DataFrame([timesr, counts]).T
    df.columns = ['Dwell Time', 'Total Dwelling Cells']
    df['Source'] = source
    sns.barplot(x='Dwell Time', y = 'Total Dwelling Cells', hue = 'Source', data = df)
    plt.title('Dwell Times')
    plt.savefig(output_location)
    plt.close()
    return

def process_dwell_times(df, output_location, exp_name,
                    n_ids=3, D=3, col_name='cell_id',
                    n_bins=15, min_count=10, max_feat=57,
                    alpha=0.05, log_trans=False):
    '''
    Calculates dwell times for each state with a given dimensional
    and binning resolution.
    '''

    df_dr = reduce_dims(df.iloc[:,n_ids:max_feat], D = D, log_trans = log_trans)
    df_id = split_col(df.iloc[:,:n_ids], col_name=col_name)
    T = int(df_id['timestep'].max())+1 # total number of timesteps
    state_locations, dim_locations = gen_state_locations(df_id, df_dr, n_bins=n_bins)
    dwell_times = gen_dwell_times(state_locations, D = D)

    sum_dwell = dwell_times.iloc[:,-T:].sum(1)
    idx = sum_dwell > min_count

    mean_dwell = calc_mean_dwell_time(dwell_times, T = T)
    plt.close('all')
    plot_dwell_v_occupancy(mean_dwell[idx], sum_dwell[idx],
                            os.path.join(output_location, exp_name + '_dwell_v_occupancy.png'))

    time_constants, scales, tc_all, scale_all = evol_dwell_time_constants(dwell_times[idx], T=T, overall_tc=True)
    plt.close('all')

    plot_tc_v_occupancy(time_constants, sum_dwell[idx],
                        os.path.join(output_location, exp_name + '_tc_v_occupancy.png'))

    # check for exponential distribution
    chi2, pval, theor = test_dwell_chisq(dwell_times, T=T)

    plt.close('all')
    plot_dwell_times(dwell_times.iloc[:,-T:].sum(0), theor, np.arange(1,T+1),
                        os.path.join(output_location, exp_name + '_dwell_times_hist.png'))

    # make a csv for output
    dwell_times_op = dwell_times[idx]
    dwell_times_op['sum_dwell'] = sum_dwell[idx]
    dwell_times_op['mean_dwell'] = mean_dwell[idx]
    dwell_times_op['time_constant'] = time_constants
    dwell_times_op['exp_chi2'] = chi2
    dwell_times_op['exp_pval'] = pval
    dwell_times_op.to_csv(os.path.join(output_location, exp_name + '_dwell_times.csv'), index=False)

    return


def process_dim_vectors(df, output_location,
                        n_ids=3, D=2,
                        col_name='cell_id', max_feat=57,
                        plot=True, size=(5,5),
                        log_trans = False, n_bins=10, w=None):
    import matplotlib.pyplot as plt
    import seaborn as sns

    df_dr = reduce_dims(df.ix[:,n_ids:max_feat], D=D, log_trans=log_trans, w=w)
    df_id = split_col(df.ix[:,:n_ids], col_name=col_name)
    state_locations, dim_locations = gen_state_locations(df_id, df_dr, n_bins=n_bins)
    dim_vectors = gen_state_vectors(dim_locations, D = D)

    if plot:
        g = sns.jointplot(x='vdim00', y='vdim01', data=dim_vectors)
        g.set_axis_labels(xlabel='dim00 vector', ylabel='dim01 vector')
        plt.savefig(output_location)

    return dim_vectors

def plot_sig_transitions(sig_p_vals, output_location=None, N=5, title=None, sigmarkers=False, cmap=None, size=(4,3)):
    '''
    Plots significantly unbalanced transitions as an N x 2
    heatmap, where N is the number of transitions plotted
    and the color intensity represents the number of cells
    undergoing a given transition.

    Parameters
    ----------
    sig_p_vals : pandas DataFrame.
        N x M dataframe, where the first n_ids columns are unique identifiers,
        with the final identifier col_name containing sample ID and timestep
        information, separated by a '-'.
        Remaining M-n_ids columns are unique features.
    output_location : string.
        output file path to save plot.
    N : integer.
        Number of the top unbalanced transitions to plot.
    sigmarkers : boolean.
        Outline significant transitions with a red box.
    cmap : seaborn compatible colormap.
        Defaults to seaborn defaults.
    size : tuple length = 2 of integers.
        Specifies x and y dimensions of the figure.

    Returns
    -------
    fig : matplotlib.pyplot figure object.
    '''
    import matplotlib.pyplot as plt
    import matplotlib.colorbar as cbar
    from matplotlib.patches import Rectangle
    import seaborn as sns

    Ntrans = sig_p_vals.iloc[:N,-2:]
    Ntrans.columns = ['t0 State', 't1 State']

    fig = plt.figure(figsize=size)
    if not cmap:
        cmap = sns.cubehelix_palette(8, reverse=False, start=2.8, rot=.1)
    g = sns.heatmap(Ntrans.astype('int32'), cbar_kws = {'label':'Cell Transition Count'}, linewidths=0, annot=True, fmt='d', cmap=cmap)

    if sigmarkers:
        ax = g.axes
        sig_idx = np.where(sig_p_vals['reject_H0'])[0]
        #rev_idx = np.tile(N-1, len(sig_idx))
        #rev_idx -= sig_idx
        rev_idx = sig_idx # diff vers of sns flip the axis!
        for y in rev_idx:
            ax.add_patch(Rectangle((0, y), 2, 1, fill = False, edgecolor='red', lw = 3))

    plt.ylabel('Transition Pair')
    if title:
        plt.title(title)

    if output_location:
        plt.savefig(output_location)

    plt.close()

    return




def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('df_location', help='input location of features with time steps')
    parser.add_argument('output_location', help='path to CSV file for output')
    parser.add_argument('--D', default=3, help='number of dimensions to consider')
    parser.add_argument('--min_count', default=6, help='min. number of transitioning cells to consider for sig. diff.')
    df = read_dataframe(df_location)
    df_dr = reduce_dims(df.ix[:,3:57], D = D) # reduce to D dims
    df_id = split_col(df.ix[:,:3], col_name='cell_id')
    state_locations = gen_state_locations(df_id, df_dr, n_bins=15)
    state_vectors = gen_state_vectors(state_locations, D = D)
    # Check if detailed balance is broken
    p_vals = pairwise_balance(state_vectors, n_bins = 15, D = D, min_count = 6)
    sig_p_vals = fdr_correction(p_vals, alpha = 0.05)
    # Write output
    sig_p_vals.to_csv(output_location, sep=',', index = False)
    print('Wrote Sig. p-values to ', output_location)
    return

if __name__ == '__main__':
    main()
