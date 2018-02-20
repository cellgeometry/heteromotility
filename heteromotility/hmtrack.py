'''
1. Establish cell paths from centroid locations and time points
2. Sanity checking of the resulting cell paths
3. Identify colliding cells, remove one or both colliders from analysis
'''
from __future__ import division, print_function
from .hmtools import dedupe, cell_ids2tracks
from .hmstats import distance, average_xy
import sys
import numpy as np

class CellPaths(object):
    '''
    Creates `cell_ids` dictionaries from a variety of inputs.

    Attributes
    ----------
    cell_ids : dict
        keyed by cell id, values are lists of coordinate tuples.
    epsilon : float.
        small float, at the magnitude of non-perturbing noise. Default = 1e-9.

    Notes
    -----
    cell_ids = {
                0 : [ (x1,y1), (x2,y2), (x3,y3), (x4,y4) ] # ea key is an object
                1 : [ (x1,y1), (x2,y2), (x3,y3), (x4,y4) ] # ea tuple is a timepoint
                }
    '''
    def __init__(self, centroid_arrays = [], cell_ids = {}, tracksX=None, tracksY=None, sanity_px = 200, interp_lim = 3):
        '''
        Creates `cell_ids` dictionaries from a variety of inputs.

        Parameters
        ----------
        centroid_arrays : list
            list of np.ndarray of centroids.
        cell_ids : dict.
            previously generated `cell_ids` dict.
        tracksX, tracksY : np.ndarray.
            N x T array of cell X, Y coordinates.
        sanity_px : int
            maximum distance for use in cell tracking with `centroid_arrays`.
        interp_lim : int
            maximum interpolation for use in cell tracking with `centroid_arrays`.
        '''

        self.epsilon = 1e-9

        if len(centroid_arrays) > 0:
            self.cell_ids, self.removed_cells, self.points_corrected = self.establish_cell_paths(centroid_arrays, sanity_px, interp_lim)
        elif len(cell_ids) > 0 and sanity_px != None:
            self.cell_ids, self.removed_cells, self.points_corrected = self.sanity_on_provided_ids(cell_ids, sanity_px)
        elif len(cell_ids) > 0 and sanity_px == None:
            self.cell_ids = cell_ids
        elif np.any(tracksX) and np.any(tracksY):
            self.cell_ids = self.tracks2cell_ids(tracksX, tracksY)
        else:
            print("CellPaths requires an array of centroids or provided cell_ids")
            sys.exit()

    def assign_cell_ids(self, centroid_arrays):
        cell_ids = {}
        i = 0
        for coor in centroid_arrays[0]:
            cell_ids[i] = []
            cell_ids[i].append(coor)
            i = i + 1
        return cell_ids

    def establish_cell_paths(self, centroid_arrays, sanity_px, interp_lim):
        cell_ids = self.assign_cell_ids(centroid_arrays)
        for array in centroid_arrays[1:]:
            for u in cell_ids:
                dists = {}
                for coor in array:
                    d = distance( coor, cell_ids[u][-1] )
                    dists[d] = coor
                closest_coor = dists[ min(dists) ]
                cell_ids[u].append(closest_coor)

        # Sanity Check
        cell_ids, remove, points_corrected = self.sanity_check(cell_ids, sanity_px, interp_lim)
        # Remove Insane
        cell_ids, removed = self.remove_cells(cell_ids, remove)

        return cell_ids, removed, points_corrected

    def sanity_on_provided_ids(self, cell_ids, sanity_px):
        cell_ids, remove, points_corrected = self.sanity_check(cell_ids, sanity_px)
        cell_ids, removed = self.remove_cells(cell_ids, remove)

        return cell_ids, removed, points_corrected

    def tracks2cell_ids(self, tracksX, tracksY):
        '''
        Converts tracks arrays to `cell_ids` dict.

        Parameters
        ----------
        tracksX : ndarray.
            N x T array of sequential X positions.
        tracksY : ndarray.
            N x T array of sequential Y positions.

        Returns
        -------
        cell_ids : dict.
            dict indexed by a unique integer for each object, containing lists
            of tuples of XY coordinates
        '''

        cell_ids = {}
        # reshape to 2D if 1D
        if len(tracksX.shape) == 1:
            tracksX = np.reshape(tracksX, (1, tracksX.shape[0]))
            tracksY = np.reshape(tracksY, (1, tracksY.shape[0]))

        N, T = tracksX.shape

        for i in range(N):
            t = []
            X = tracksX[i,:]
            Y = tracksY[i,:]
            for j in range(T):
                t.append((X[j], Y[j]))
            cell_ids[i] = t

        return cell_ids

    def sanity_check(self, cell_ids, sanity_px, interp_lim = 3):
        points_corrected = []
        remove = []
        for u in cell_ids:
            cell = cell_ids[u]

            for coor in cell[1:-4]: # exclude the first & last 4 coordinate
                previous_coor = cell[ ( cell.index(coor) - 1 ) ]
                d = distance( coor, previous_coor )
                if d > sanity_px:
                    i = 1
                    fix = False
                    next_coor = cell[ (cell.index(coor) + i) ]
                    while i < (interp_lim + 1):
                        d2 = distance( previous_coor, next_coor )
                        if d2 < sanity_px:
                            fixed_coor = average_xy( previous_coor, next_coor )
                            cell_ids[u][ cell.index(coor) ] = fixed_coor
                            correction = [ coor, fixed_coor, previous_coor, next_coor ]
                            points_corrected.append(correction)
                            fix = True
                            i = interp_lim + 1
                        else:
                            i += 1
                    if fix == False:
                        remove.append(u)
                    else:
                        pass
                else:
                    pass

        # Remove duplicate removal marks from a single cell
        remove = dedupe(remove)

        return cell_ids, remove, points_corrected

    def remove_cells(self, cell_ids, to_remove):
        removed = {}
        for u in to_remove:
            removed[u] = cell_ids.pop(u, None)

        return cell_ids, removed

    def check_invariant_tracks(self):
        '''
        Checks for non-zero movement in both `X` and `Y` dimensions.
        '''
        invariant = np.zeros((len(self.cell_ids), 2))

        i = 0
        for u in self.cell_ids:
            path = self.cell_ids[u]
            XY = np.array(path)
            x_unique = np.unique(XY[:,0])
            y_unique = np.unique(XY[:,1])
            invariant[i,:] = (len(x_unique) <= 1, len(y_unique) <= 1)
            i += 1

        self.invariant_tracks = invariant
        return

    def fuzz_invariant_tracks(self):
        '''
        Applies random epsilon noise to prevent invariant tracks
        if an invariant track has been detected
        '''
        if self.invariant_tracks.sum() > 0:

            tracksX, tracksY = cell_ids2tracks(self.cell_ids)
            tracksX += (2+np.random.randn(*tracksX.shape))*self.epsilon
            tracksY += (2+np.random.randn(*tracksY.shape))*self.epsilon
            tracksX[tracksX<0] = 0
            tracksY[tracksY<0] = 0 # in case fuzz made negative vals
            self.cell_ids = self.tracks2cell_ids(tracksX, tracksY)

        return

def check_remaining_cells( cell_ids, verbose=False ):

    if len(cell_ids) > 0:
        if verbose:
            print(len(cell_ids), ' cells remain after sanity checking \n')
    else:
        if verbose:
            print("No cells remain after sanity checking \n")
        sys.exit()
