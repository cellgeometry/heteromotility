from .hmstats import GeneralFeatures, MSDFeatures, RWFeatures
from .hmtrack import CellPaths
from .hmio import *
from .hmtools import *

class DetailedBalance(object):
    '''
    Functions to split cell path's into discrete segments and calculate
    motility statistics on these sub-segments,
    such that 'motility state' transitions may be assessed within a single cell
    and the relevant timescale for 'motility state' defintion may be defined

    '''
    def __init__(self, cell_ids, min_split = 20, tau_max=None):
        #self.multi_split = self.multi_tau_split(cell_ids, tau_min = min_split, tau_max = tau_max)
        self.min_split = min_split
        self.tau_max = tau_max
        self.cell_ids = cell_ids
        pass

    def split_ids(self, cell_ids, tau):
        '''
        Splits cell paths into multiple discrete segments of length tau
        Returns a dict of lists where each list is a subpath for easy stat
        calculations with existing Heteromotility infrastructure

        Parameters
        ----------------------
        cell_ids : dict of lists, ea. list containing tuples of XY coordinates
        tau : time lag, size of each split path

        Returns
        ----------------------
        split_ids :
        dict of lists, keyed by "cell_id-subpath_id"
        '''

        split_ids = {}
        n = len(cell_ids[ list(cell_ids)[0] ])
        num_subseries = int(n/tau)
        for u in cell_ids:
            path = cell_ids[u]
            for i in range(0, num_subseries):
                sub = path[ i*tau : (i*tau) + tau ]
                split_ids[str(u) + '-' + str(i)] = sub

        return split_ids

    def sliding_windows(self, cell_ids, tau=20, stride=1, max_windows=None):
        '''
        Generates multiple tracks using a sliding window of length `tau`,
        shifted by `stride` for each frame until the track is finished or
        `max_windows` is reached.

        Parameters
        ----------
        cell_ids : dict.
            dict of lists, ea. list containing tuples of XY coordinates.
        tau : integer.
            size of the sliding window in time step units.
        stride : integer.
            stride of the sliding window in time step units.
        max_windows : integer.
            maximum number of sliding windows to extract.
            if `None`, extracts the maximum number of windows across the track.

        Returns
        -------
        sliding_ids : dict.
            dict of lists, keyed by `"cell_id-subpath_id"` where `subpath_id` is
            equal to the window location, indexed from `0`.
        '''

        sliding_ids = {}
        T = len(cell_ids[ list(cell_ids)[0] ]) # track lengths
        if not max_windows:
            # output = [ floor(input - kernel) / stride ] + 1
            num_windows = [ np.floor(T - tau) / stride ] + 1
        else:
            num_windows = max_windows

        for u in cell_ids:
            start = 0
            for i in range(0, num_windows):
                path = cell_ids[u]
                window = path[start : start + tau]
                sliding_ids[str(u) + '-' + str(i)] = window
                start += stride

        return sliding_ids

    def multi_tau_split(self, cell_ids, tau_min = 20, tau_max = None):
        '''
        Makes dict keyed by tau with split_ids for a range of time lags

        Parameters
        ----------------------
        cell_ids : dict of lists, ea. list containing tuples of XY coordinates
        tau_min : minimum time lag to consider

        n.b. max_tau is internally set at 0.5*len(time_series)

        Returns
        ----------------------
        '''
        if tau_max == None:
            tau_max = len( cell_ids[ list(cell_ids)[0] ] )
            tau_max = int(tau_max / 2)

        tau_range = range( tau_min, tau_max + 1)
        multi_split = {}
        for tau in tau_range:
            multi_split[tau] = self.split_ids(cell_ids, tau)

        return multi_split

    def split_id_features(self, multi_split, output_dir, output_suffix = False):
        '''
        Calculates motility features for subsegments of cell
        paths. Exports data as 'motility_statistics_split_$TAU.csv'
        where $TAU is the size of the segments being considered

        Parameters:
        multi_split : dict keyed by tau, containing cell_id style dicts
        output_dir : directory for output files
        '''
        for tau in multi_split:
            print("Detailed Balance is calculating stats for paths of length :", tau)
            split_ids = multi_split[tau]
            scp = CellPaths(cell_ids = split_ids, sanity_px = None)
            sgf = GeneralFeatures(cell_ids = split_ids, move_thresh = 10)
            smsdf = MSDFeatures(cell_ids = split_ids)
            srwf = RWFeatures(cell_ids = split_ids, gf= sgf)

            ind_outputs = single_outputs_list(scp.cell_ids, sgf, srwf, smsdf, output_dir, suffix=output_suffix)
            merged_list = make_merged_list(ind_outputs, sgf, srwf)
            if type(output_suffix) != str:
                output_name = 'motility_statistics_split_' + str(tau) + '.csv'
            else:
                output_name = 'motility_statistics_split_' + str(tau) + '_' + output_suffix + '.csv'
            write_motility_stats(output_dir, output_name, sgf, srwf, merged_list)


        return merged_list
