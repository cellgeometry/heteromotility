from hmtrack import *

'''
#---------------------------
# MODULE CONTENTS
#---------------------------
Functions to/for:

Test efficacy of heteromotility package

'''
#------------------------------
# Tests
#------------------------------
# Tests cell_id generation

# Tests feature calculation functions

# Tests sanity_check corrections

def test_sanity():
	sanity_px = 5

	test_data = { 0:[ (1, 2), (15, 20), (1, 3), (2, 3), (2,3), (2,3), (2,4), (2,3), (2,4)],
  			      1:[ (1, 1), (1, 3), (1, 5), (0, 5), (1,5), (1,5), (1,5), (1,6), (2,6)],
				  2:[ (4, 3), (20, 43), (4, 4),  (5, 4), (5,5), (5,5), (6,5), (5,6), (5,5)]}

	proper_ids = {0: [(1, 2), (1.0, 2.5), (1, 3), (2, 3), (2, 3), (2, 3), (2,4), (2,3), (2,4)],
				  1: [(1, 1), (1, 3), (1, 5), (0, 5), (1, 5), (1, 5), (1,5), (1,6), (2,6)],
				  2: [(4, 3), (4.0, 3.5), (4, 4), (5, 4), (5, 5), (5, 5), (6,5), (5,6), (5,5)]}
	proper_removed = {}
	proper_changes =  [ [(15, 20), (1.0, 2.5), (1, 2), (1, 3)],
					    [(20, 43), (4.0, 3.5), (4, 3), (4, 4)] ]

	cp = CellPaths(cell_ids = test_data, sanity_px = 5)
	test_ids = cp.cell_ids
	test_removed = cp.removed_cells
	test_changes = cp.points_corrected

	if proper_ids == test_ids and proper_removed == test_removed and proper_changes == test_changes:
		return True
	else:
		return False

# Test removal of lost cells

def test_removal():

	sanity_px = 5
	test_data = { 0:[ (1, 2), (1 ,3), (15, 20), (10, 30), (25, 35), (12,32), (22,3), (20, 20), (1,5), (1,4)],
  			      1:[ (1, 1), (1, 3), (1, 5), (0, 5), (1,5), (1,5), (1,6), (1,7)],
				  2:[ (4, 3), (4, 3), (20, 43), (14, 14),  (25, 4), (45,5), (15,15), (20, 20), (4,4), (4,3)]}

	proper_ids = { 1:[ (1, 1), (1, 3), (1, 5), (0, 5), (1,5), (1,5), (1,6), (1,7)] }

	proper_removed = [0,2]
	proper_changes = [[(25, 35), (11.0, 31.0), (10, 30), (12, 32)]]

	cp = CellPaths(cell_ids = test_data, sanity_px = 5)
	test_ids = cp.cell_ids
	test_removed = list(cp.removed_cells)
	test_changes = cp.points_corrected

	if proper_ids == test_ids and proper_removed == test_removed and proper_changes == test_changes:
		result =  True
	else:
		result = False

	return result #,test_ids, test_removed, test_changes

# Unit test for finding colliders
def test_colliders():

	test_data = {
	0:[(0,0),(1,3),(1,5),(2,7),(3,5),(6,6),(5,5),(6,6),(5,5),(6,6),(5,5),(7,7)],
	1:[(5,4),(3,2),(3,3),(5,5),(6,6),(3,5)],
	2:[(3,3),(2,3),(1,5),(2,7),(3,4),(3,4)],
	3:[(2,4),(5,4),(3,4),(5,5),(3,2),(6,7)],
	4:[(3,3),(7,8),(8,9),(5,5),(2,1),(9,9)]
	}

	colliders, colliding_t, prop_colliders = find_colliders(test_data)

	correct_colliders = [2,2]
	correct_coliding_t = [2,3]
	correct_prop_colliders = ( 2/3 )

	correct = [correct_colliders, correct_coliding_t, correct_prop_colliders]
	test = [colliders, colliding_t, prop_colliders]

	if test == correct:
		return True
	else:
		return False
