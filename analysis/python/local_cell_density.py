'''
Local cell density

Determine the local cell density within a given imaging frame
'''

import numpy as np
import os
import pandas
import glob
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix

def density_at_t(pos):
    '''
    Calculates sum of distances to other cells.

    Parameters
    ----------
    pos : ndarray. N x 2, [x, y] coordinates for each row.

    Returns
    -------
    distances : ndarray. N x 1.
        sum of (neighboring cells / distance**2).
    '''
    dm = distance_matrix(pos, pos)
    dm[dm==0] = np.inf # set zeros to inf for inversion
    inv = dm**-2
    s = inv.sum(1)
    return s

def density_over_time(x, y):
    '''
    Calculates distance to neighboring cells for each time point in a series.

    Parameters
    ----------
    x : ndarray. N x T.
        x coordinates.
    y : ndarray. N x T.
        y coordinates.

    Returns
    -------
    dist_timeseries : ndarray. N x T.
        sums of distances to neighboring cells.
    '''
    assert x.shape == y.shape, 'x and y are not the same shape'

    if len(x.shape) == 1:
        x = np.reshape(x, (1, x.shape[0]))
        y = np.reshape(y, (1, y.shape[0]))

    T = x.shape[1]
    dist_timeseries = np.zeros(x.shape)
    for t in range(T):
        pos = np.stack([x[:,t], y[:,t]]).T
        dist_timeseries[:,t] = density_at_t(pos)

    return dist_timeseries

def mean_density_values(xy, in_dir, out_dir):
    '''
    Saves mean density over time as a CSV.
    '''

    try:
        x = np.loadtxt(os.path.join(in_dir, 'tracksX_xy' + str(xy).zfill(2) + '.csv'), delimiter = ',')
        y = np.loadtxt(os.path.join(in_dir, 'tracksY_xy' + str(xy).zfill(2) + '.csv'), delimiter = ',')
    except:
        print('File not found: ', os.path.join(in_dir, 'tracksX_xy' + str(xy).zfill(2) + '.csv'))
        return

    if np.sum(x) == 0:
        print('File empty: ', os.path.join(in_dir, 'tracksX_xy' + str(xy).zfill(2) + '.csv'))
        return

    dens_timeseries = density_over_time(x,y)
    means = dens_timeseries.mean(1)

    np.savetxt(os.path.join(out_dir, 'lcdensity_xy' + str(xy).zfill(2) + '.csv'), means, delimiter = ',')
    return

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('in_dir', help='input directory containing tracksX_xy/tracksY_xy')
    parser.add_argument('--out_dir', type=str, default='', help='output directory. Default = input dir.')
    parser.add_argument('--max_xy', type=int, default=25, help='maximum number of xy positions. Default = 25.')
    args = parser.parse_args()

    in_dir = args.in_dir
    if args.out_dir == '':
        out_dir = args.in_dir
    else:
        out_dir = args.out_dir
    max_xy = args.max_xy

    for xy in range(max_xy):
        print('Local density calculation ', in_dir, xy+1)
        mean_density_values(xy+1, in_dir, out_dir)
    return

if __name__ == '__main__':
    main()
