#------------------------------
# PLOT CELL MOTILITY PATHS
#------------------------------

def plot_paths(cell_ids):
	for u in cell_ids:
		cell = cell_ids[u]
		x = []
		y = []
		for coor in cell:
			x.append( float(coor[0]) )
			y.append( float(coor[1]) )

	# Draws a scatter plot of x, y values for each path
	# Connects plots with lines
		plt.scatter(x, y)
		plt.plot(x, y)



#------------------------------
# FIT SPLINES TO CELL PATHS
#------------------------------

def fit_splines(cell_ids):
	splines = [] # list of np.arrays in tuple pairs (x_array, y_array)
	for u in cell_ids:
		cell = cell_ids[u]
		x = []
		y = []
		for coor in cell:
			x.append( float(coor[0]) )
			y.append( float(coor[1]) )

		# Calculate a range spline
		x = np.asarray(x)
		y = np.asarray(y)
		t = np.arange(x.shape[0], dtype=float)
		t /= t[-1]
		nt = np.linspace(0, 1, 100)
		x1 = interpolate.spline( t, x, nt )
		y1 = interpolate.spline( t, y, nt )

		# Calculate a distance spline
		t = np.zeros(x.shape)
		t[1:] = np.sqrt( (x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2 )
		t = np.cumsum(t)
		t /= t[-1]
		x2 = interpolate.spline(t, x, nt)
		y2 = interpolate.spline(t, y, nt)

		splines.append( (x1, y1 , x2, y2) )


	return splines, x1, y1, x2, y2
