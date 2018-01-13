from __future__ import division, print_function
from .hmtools import dedupe, cell_ids2tracks
from .hmstats import distance, average_xy
import sys
import numpy as np

'''
#---------------------------
# MODULE CONTENTS
#---------------------------
Classes to/for:
1. Establish cell paths from centroid locations and time points
2. Sanity checking of the resulting cell paths
3. Identify colliding cells, remove one or both colliders from analysis

'''

#---------------------------
# ESTABLISH CELL PATHS
#---------------------------
# WIP
class CellPaths:
    '''
    CellPaths accepts either a list of lists of tuples containing centroid
    XY coordinates for each object in a frame (centroid_arrays)
    OR
    a dict keyed by unique cell identifiers
    (i.e. 0, 1, 2) with values as a list of tuples containg a single cell's XY
    coordinates, each tuple representing a given timepoint

    If an array of centroids is provided, CellPaths uses Heteromotility's
    built in tracking algorithm
    A dict of lists two ndarrays of X and Y coordinates would be provided
    if another tracking algorithm is desired

    centroid_arrays = [
                        [(x,y), (x,y), (x,y), ...] # ea. sublist is a timepoint
                        [(x,y), (x,y), (x,y), # each tuple is one object's XY ]
                        ]

    cell_ids = {
                0 : [ (x1,y1), (x2,y2), (x3,y3), (x4,y4) ] # ea key is an object
                1 : [ (x1,y1), (x2,y2), (x3,y3), (x4,y4) ] # ea tuple is a timepoint
                }
    '''
    def __init__(self, centroid_arrays = [], cell_ids = {}, tracksX=None, tracksY=None, sanity_px = 200, interp_lim = 3):

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

        # self.check_invariant_tracks()
        # self.fuzz_invariant_tracks()

    # assign each object in t0 centroids array a key in cell_ids dict
    # keys are just iterated integers, could do something fancier
    # if you need to for any reason
    def assign_cell_ids(self, centroid_arrays):
        cell_ids = {}
        i = 0
        for coor in centroid_arrays[0]:
            cell_ids[i] = []
            cell_ids[i].append(coor)
            i = i + 1
        return cell_ids

    # for each unique cell in cell_ids, find the closest coordinate to
    # the previous value in the next centroid array and append it
    # NOTE:
    # This assumes fairly limited cell motility b/w frames
    # and very little cell::cell contact

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

    #---------------------------
    # SANITY TESTING OF PATHS
    #---------------------------

    # Find points that move too quickly
    # (likely not segmented in one frame)
    # Correct for this by interpolating from the previous sensible point
    # and the next sensible point
    # Outputs a corrected cell_ids dict, a list of cells to remove,
    # and a list of lists with the points that were corrected

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


#------------------------------
# CHECK FOR REMAINING CELLS
#------------------------------

# Checks to see if any cells remain after removing cells
# that fail the sanity check
# If no cells are left, exits the script gracefully
# without pickling cell_ids or removed_ids

def check_remaining_cells( cell_ids, verbose=False ):

    if len(cell_ids) > 0:
        if verbose:
            print(len(cell_ids), ' cells remain after sanity checking \n')
    else:
        if verbose:
            print("No cells remain after sanity checking \n")
        sys.exit()

#------------------------------
# CHECK TWO CELLS COLLIDING
#------------------------------


# Detects collisions of cells, and records them in a dict collider_log
# Stores all seen positions and cooresponding cells in seen = {}

# seen = { cell1 : coor, cell2 : coor, ... }

# When a collision is found, adds the colliding cells and coor to
# collider log

# collider_log = { t1: [ [cell1, cell2, coor], [cell3, cell4, coor]
#                   t2: ...                                           }

def find_colliders(cell_ids):

    collider_log = {}
    colliding_cells_log = {}
    seen_log = {}
    t = 0
    first_key = list(cell_ids)[0]
    while t < len( cell_ids[ first_key ] ):
        seen = {}
        colliders = []
        colliding_cells = []

        for u in cell_ids:
            coor = cell_ids[u][t]
            if coor in seen.values():

                for prev_cell, prev_coor in seen.iteritems():
                    if coor == prev_coor:
                        colliders.append( [ (prev_cell, u) , coor] )
                        colliding_cells.append(u)
                        colliding_cells.append(prev_cell)
                    else:
                        pass
                seen[u] = coor
            else:
                seen[u] = coor

        seen_log[t] = seen
        collider_log[t] = colliders
        colliding_cells_log[t] = dedupe(colliding_cells)
        prop_colliders = float(len(dedupe(colliding_cells))) / float(len(cell_ids))

        t += 1

    return collider_log, prop_colliders

# Takes collider_log, checks to see if collisions
# are longer than a stringencey coeff.
# If not, marks cells as persistent colliders

def find_persistent_colliders(collider_log, stringency_coeff):

    t = 0
    persistent_colliders = []
    temp_colliders = []
    while t < (len(collider_log) - 30):
        timepoint = collider_log[t]
        for group in timepoint:
            collide_pair = group[0]
            future_pairs = []
            for i in range(1,31):
                future = collider_log[ t + i ]
                for future_group in future:
                    future_pairs.append( future_group[0] )

            if future_pairs.count( collide_pair ) > stringency_coeff:
                persistent_colliders.append( collide_pair )
            else:
                temp_colliders.append( collide_pair )
        t += 1

    persistent_colliders = dedupe(persistent_colliders)
    temp_colliders = dedupe(temp_colliders)

    return persistent_colliders, temp_colliders

# Takes persistent colliders & uses removal function
# from sanity testing to remove from cell_ids

def remove_both_colliders(persistent_colliders, cell_ids):

    to_remove = []
    for collide_pair in persistent_colliders:
        to_remove.append(collide_pair[0])
        to_remove.append(collide_pair[1])

    to_remove = dedupe(to_remove)
    to_remove.sort()

    cell_ids, collider_ids = remove_cells(cell_ids, to_remove)
    return cell_ids, collider_ids

def remove_one_collider(persistent_colliders, cell_ids):

    to_remove = []
    for collide_pair in persistent_colliders:
        to_remove.append(collide_pair[0])

    to_remove = dedupe(to_remove)
    to_remove.sort()

    cell_ids, collider_ids = remove_cells(cell_ids, to_remove)
    return cell_ids, collider_ids

def check_colliding_cells( cell_ids ):

    if len(cell_ids) > 0:
        print(len(cell_ids), ' cells remain after collider checking \n')
    else:
        print("No cells remain after collider checking \n")
        sys.exit()
