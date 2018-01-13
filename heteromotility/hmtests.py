'''
Heteromotilit test suite
'''
import numpy as np
from numpy.testing import assert_almost_equal
import pandas as pd

from .hmtrack import CellPaths
from .hmstats import GeneralFeatures, MSDFeatures, RWFeatures
from .hmtools import make_merged_list, single_outputs_list

class TestInvariantTrackFuzzing(object):
    '''
    Tests fuzzing of tracks which have an invariant dimension
    '''

    def test_check_invariant(self):

        x = np.zeros([10,50])
        y = np.zeros([10,50])
        fuzz = np.random.random(x.shape)
        fuzz[9,:] = 0
        x += fuzz
        fuzz = np.random.random(x.shape)
        fuzz[7,:] = 0
        y += fuzz

        cp = CellPaths(tracksX=x, tracksY=y)
        idx = np.where(cp.invariant_tracks.sum(1)>0)[0]
        assert (7 in idx)
        assert (9 in idx)

    def test_fuzzed_stats(self):

        x = np.zeros((3,50))
        x[0,:] = np.linspace(1,10,50)
        x[1,:] = np.linspace(17,25,50)
        x[2,:] = np.linspace(9,33,50)

        y = np.zeros((3,50))
        y[0,:] = np.linspace(0,5,50)
        y[2,:] = np.linspace(11,21,50)

        cp0 = CellPaths(tracksX=x, tracksY=y)
        cp1 = CellPaths(tracksX=x, tracksY=y)

        gf0 = GeneralFeatures(cp0.cell_ids, move_thresh = 0.1)
        msdf0 = MSDFeatures(cp0.cell_ids)
        rwf0 = RWFeatures(cp0.cell_ids, gf0)
        ind_outputs0 = single_outputs_list(cp0.cell_ids, gf0, rwf0, msdf0, '.', suffix='')
        merged_list0 = make_merged_list(ind_outputs0, gf0, rwf0)

        gf1 = GeneralFeatures(cp1.cell_ids, move_thresh = 0.1)
        msdf1 = MSDFeatures(cp1.cell_ids)
        rwf1 = RWFeatures(cp1.cell_ids, gf1)
        ind_outputs1 = single_outputs_list(cp1.cell_ids, gf1, rwf1, msdf1, '.', suffix='')
        merged_list1 = make_merged_list(ind_outputs1, gf1, rwf1)

        df0 = pd.DataFrame(merged_list0)
        df1 = pd.DataFrame(merged_list1)
        # delete undefined features for invariant tracks
        # linearity, spearman, etc.
        drop_idx = [0, 1, 4, 5, 15] + list(range(47, df0.shape[1]))
        df0 = df0.drop(drop_idx, axis=1)
        df1 = df1.drop(drop_idx, axis=1)

        diff = df0 - df1
