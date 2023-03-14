'''
1. Calculate general motility features (average speed, distance, ...)
2. Calculate features related to mean squared displacement
3. Calculate features related to random walk similarity
'''

from __future__ import division, print_function
import numpy as np
import math
from scipy import stats
from scipy import interpolate
from statsmodels.tsa.stattools import acf, pacf

def average_xy( coor1, coor2 ):
    '''Finds the average value of two XY coordinates'''
    x1 = float(coor1[0])
    x2 = float(coor2[0])
    y1 = float(coor1[1])
    y2 = float(coor2[1])

    x3 = (x1 + x2) * 0.50
    y3 = (y1 + y2) * 0.50
    return (x3, y3)

def distance( coor1, coor2 ):
    '''Euclidean distance'''
    x1, y1 = coor1
    x2, y2 = coor2
    d = math.sqrt( (x2-x1)**2.00 + (y2-y1)**2.00)
    return d

class GeneralFeatures(object):
    '''
    Calculates features of speed, distance, path shape, and turning behavior.

    Attributes
    ----------
    cell_ids : dict
        keyed by cell ids, values are lists of coordinate tuples.
    moving_tau : int
        threshold for movement
    tau_range : iterable, ints
        range of window sizes to use to establish cell direction by regression.
    interval_range : iterable, ints
        range of time intervals to use between a cell's current and 'next'
        direction when calculating turning features.
    total_distance : dict
        total distance each cell travels.
        keyed by cell id, values are floats.
    net_distance : dict
        net distance each cell traveled, dist(final_pos - init_pos)
        keyed by cell id, values are floats.
    linearity : dict
        Pearson's r**2 linearity metric on cell paths.
        keyed by cell id, values are floats.
    spearmanrsq : dict
        Spearman's rho**2 monotonicity metric on cell paths.
        keyed by cell id, values are floats.
    progressivity : dict
        net_distance/total_distance, metric of directional persistence.
        keyed by cell id, values are floats.
    min_speed : dict
        minimum cell speed.
        keyed by cell id, values are floats.
    max_speed : dict
        maximum cell speed.
        keyed by cell id, values are floats.
    avg_moving_speed : dict
        average moving speed above a certain threshold speed,
        defined as the movement cutoff.
        useful to deconfound the cell's average speed while moving from the
        dwell time in an immotile state.
        keyed by threshold, values are dicts keyed by cell id, valued with floats.
    time_moving : dict
        proportion of time spent moving above a certain threshold speed,
        defined as the movement cutoff.
        keyed by threshold, values are dicts keyed by cell id, valued with floats.
    turn_stats : dict
        proportion of the time a cell turns to the right of its past direction.
        useful to determine directional bias.
        keyed by `tau` from `tau_range`, values are dicts keyed by `interval` from `interval_range`.
        final dicts have float values, [0, 1].
    theta_stats : dict
        average turning angle magnitude.
        keyed by `tau` from `tau_range`, values are dicts keyed by `interval` from `interval_range`.
        final dicts have float values, [0, 1].
    '''
    def __init__(self, cell_ids, move_thresh = 5):
        '''

        Parameters
        ----------
        cell_ids : dict
            keyed by cell ids, values are lists of coordinate tuples.
        move_thresh : int
            time interval to use when calculating `avg_moving_speed`, `time_moving`
        '''
        self.cell_ids = cell_ids
        if len(cell_ids[list(cell_ids)[0]]) > 19:
            self.moving_tau = move_thresh
            self.tau_range = range(9,12)
            self.interval_range = range(5,7)
        elif 9 < len(cell_ids[list(cell_ids)[0]]) <= 19:
            self.moving_tau = move_thresh
            self.tau_range = range(6,7)
            self.interval_range = range(3,4)
        else:
            print("Time series too small to calculate turning features")
        self.total_distance, self.avg_speed, self.d = self.calc_total_distance(self.cell_ids)
        self.net_distance = self.calc_net_distance(self.cell_ids)
        self.linearity = self.calc_linearity(self.cell_ids)
        self.spearmanrsq = self.calc_spearman(self.cell_ids)
        self.progressivity = self.calc_progressivity(self.cell_ids, net_distance=self.net_distance, total_distance=self.total_distance)
        self.min_speed, self.max_speed = self.calc_minmax_speed(self.cell_ids)
        self.avg_moving_speed, self.time_moving = self.calc_moving_variable_threshold(cell_ids, thresholds = range(1,11), tau = self.moving_tau)

        self.turn_stats, self.theta_stats = self.all_turn_stats(cell_ids, tau_range = self.tau_range, interval_range = self.interval_range)

    def calc_total_distance(self, cell_ids):
        total_distance = {}
        avg_speed = {}
        distance_dist = {}
        for u in cell_ids:
            cell = cell_ids[u]
            start_point = cell[0]
            traveled = []
            for coor in cell[1:]:
                d = distance( start_point, coor )
                traveled.append(d)
                start_point = coor
            total_distance[u] = sum(traveled)
            avg_speed[u] = np.mean(traveled)
            distance_dist[u] = traveled

        return total_distance, avg_speed, distance_dist

    # Calculate net distance traveled and save each cell's stat
    # as a float in a list net_distance
    # Simple distance b/w start and end points
    def calc_net_distance(self, cell_ids):
        net_distance = {}
        for u in cell_ids:
            cell = cell_ids[u]
            start_point = cell[0]
            end_point = cell[-1]
            d = distance( start_point, end_point)
            net_distance[u] = d
        return net_distance

    # Calculate linearity of the traveled path as the R^2 value
    # of a linear regression of all path points
    def calc_linearity(self, cell_ids):
        linearity = {}
        for u in cell_ids:
            cell = cell_ids[u]
            x = []
            y = []
            for coor in cell:
                x.append( coor[0] )
                y.append( coor[1] )
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            r_squared = r_value**2
            linearity[u] = r_squared
        return linearity

    # Calculate Spearman's rho^2 for all XY points in a cell path
    def calc_spearman(self, cell_ids):
        spearmanrsq = {}
        for u in cell_ids:
            xy = np.array( cell_ids[u] )
            if len(np.unique(xy[:,0])) == 1 or len(np.unique(xy[:,1])) == 1:
                # fuzz to avoid bounds issue
                print('Tracks for cell_id %s were invariant in at least one dimension.' % u)
                print('Spearman coefficient is undefined, and should not be used for analysis.')
                rho = 0.
            else:
                rho, pval = stats.spearmanr( xy )
            spearmanrsq[u] = rho**2
        return spearmanrsq

    def calc_progressivity(self, cell_ids, net_distance, total_distance):
        progressivity = {}
        for u in cell_ids:
            progressivity[u] = float(net_distance[u]) / float(total_distance[u])

        return progressivity

    def calc_minmax_speed(self, cell_ids):
        '''
        Calculate max/min speeds as the distance traveled
        across 5 frames and append as a float to list max_speed/min_speed
        Speeds in the unit of pixels per hour
        '''
        max_speed = {}
        min_speed = {}
        for u in cell_ids:
            cell = cell_ids[u]
            speeds = []
            i = 0
            while (i + 4) < len(cell):
                d = distance( cell[i], cell[i+4] )
                speeds.append( d / 5 )
                i += 1
            max_speed[u] = max(speeds)
            min_speed[u] = min(speeds)
        return min_speed, max_speed

    def calc_moving_stats(self, cell_ids, move_speed, tau):
        '''
        Calculates the proportion of time a cell spends moving_speeds,
        min, max, and average speeds while moving

        Parameters
        ------------
        cell_ids : dict keyed by cell id, containing lists of tupled XY coors
        move_speed : minimum speed to be considered "moving"
        tau : time lag to use when determining movement speed
        '''
        avg_moving_speed = {}
        time_moving = {}
        for u in cell_ids:
            cell = cell_ids[u]
            speeds = []
            i = 0
            while (i + tau) < len(cell):
                d = distance( cell[i], cell[i+tau] )
                speeds.append( d / tau )
                i += 1
            moving_speeds = []
            for val in speeds:
                if val > move_speed:
                    moving_speeds.append( val )
                else:
                    pass

            time_moving[u] = float(len(moving_speeds)) / float(len(speeds))

            if moving_speeds == []:
                moving_speeds.append(0)
            else:
                pass
            avg_moving_speed[u] = np.mean( moving_speeds )

        return avg_moving_speed, time_moving

    def calc_moving_variable_threshold(self, cell_ids, thresholds = range(1,11), tau = 5):

        # dicts keyed by moving speed threshold used for calc
        # each key corresponds to a dict keyed by cell_id, with vals for moving speed
        # and time moving for the given threshold in the top level dict
        avg_moving_speed = {}
        time_moving = {}

        for thresh in thresholds:
            avg_moving_speed[thresh], time_moving[thresh] = self.calc_moving_stats( cell_ids, thresh, tau )

        return avg_moving_speed, time_moving

    def turn_direction(self, cell_path, tau = 10, interval = 5):
        '''
        Returns the proportion of time a cell turns left as float 0.0-1.0.

        Parameters
        -----------
        cell_path : list
            timeseries of tuples, each tuple holding an XY position
            i.e. cell_path = [(x1,y1), (x2,y2)...(xn,yn)]
        tau : int
            desired time lag between a point of
            interest and a point in the distance (p_n+tau) to determine
            a cell's turning behavior. Must be > `interval`.
        interval  : int
            number of points ahead and behind the point of interest to consider
            when find a regression line to represent the cell's direction at p_n.

        Returns
        -------
        turns  : list
            binaries reflecting turns left (0) or right (1)
        thetas : list
            angles of each turn in radians

        Notes
        -----
        (1) Estabish direction of object at p_n:
        Given a time series XY of x,y coordinates, where p_t denotes a point
        at time t, take a given point p_n
        determine the 'direction' of motion at p_n by plotting a linear reg.
        on points [p_n-i, p_n+i] for some interval value i
        This linear regression function is R(x) = slope*x + b

        (2) Determine if p_n+tau is left or right
        For a point p_n+tau, for a variable time lag tau, determine if p_n+tau
        lies left or right of R
            p_n+tau,Y > R(p_n+tau) and initial direction = right : left turn
            p_n+tau,Y < R(p_n+tau) and initial direction = right : right turn
            ...

        (3) Calculate magnitude of turn theta
        Calculate the angle of the turn theta as the angle between
            R(x) from (1)
            R2, line connecting p_n to p_n+tau

            v1 = (1, slope of R(x))
            v2 = (1, slope of R2)

            theta = arccos( v1 dot v2 / ||v1|| * ||v2||)

        N.B. If a cell moves in a perfectly straight line, such that point of
        interest p_n has x coordinate == p_n+tau, the function skips this p_n
        i.e. Skips if the cell has moved in a perfectly linear line, or circled
        back on itself
        '''
        if tau < interval:
            print('tau must be > interval')
            return None
        else:
            pass

        XY = np.asarray(cell_path)

        turns = []
        thetas = []
        # do every 3rd point to cut runtime
        for n in range( interval, len(XY)-tau ):
            if len(XY) > 30 and n % 3 != 0:
                continue
            p_n = XY[n]
            p_tau = XY[n+tau]
            if p_n[0] == p_tau[0]:
                continue
            direction = p_tau[0] - p_n[0]
            regress_p = XY[n-interval:n+interval]
            # skip p_n if all x coors in interval are equal
            if all(regress_p[:,0] == regress_p[0,0]):
                continue
            R_slope, b, rval, pval, stderr = stats.linregress(regress_p[:,0], regress_p[:,1])
            # skip if cell moves entirely on a vertical line for given interval period
            if np.isnan(R_slope):
                continue
            R_tau = R_slope*(p_tau[0]) + b
            # direction > 0 = +x movement, < 0 = -x
            # turn = 0 : left
            # turn = 1 : right
            if direction > 0:
                if R_tau < p_tau[1]:
                    turn = 0
                else:
                    turn = 1
            else:
                if R_tau < p_tau[1]:
                    turn = 1
                else:
                    turn = 0
            turns.append(turn)

            # find line R2 connecting p_n and p_tau
            R2_slope = (p_tau[1] - p_n[1]) / (p_tau[0] - p_n[0])
            R2_intercept = p_n[1] - R2_slope*p_n[0]
            # create vectors v1 : R, v2 : R2
            v1 = np.array( [1, R_slope] )
            v2 = np.array( [1, R2_slope] )
            # cos(theta) = v1 dot v2 / ||v1|| * ||v2||
            if (np.linalg.norm(v1)*np.linalg.norm(v2)) == 0:
                print("error (np.linalg.norm(v1)*np.linalg.norm(v2)) == 0")
                print(n, p_n, p_tau)
                continue
            costheta = np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2))
            #print "v1 = ", v1, "| v2 = ", v2, "| costheta = ", costheta
            if costheta < -1.0 or costheta > 1.0:
                #print "Error calculating turn magnitude"
                #print "costheta = ", costheta
                #print "v1, v2 = ", v1, v2
                continue
            theta = math.acos(costheta)
            thetas.append(theta)

            # If no turns are detected, add neutral values to avoid np.mean Error
            # 0 turn mag, neutral turn direction

        if len(thetas) == 0:
            thetas.append(0)
            turns.append(0.5)
            return turns, thetas
        else:
            return turns, thetas

    def all_turn_stats(self, cell_ids,
                        tau_range = range(9,12), interval_range = range(5,7)):
        '''
        Parameters
        ----------
        cell_ids : dict
            keyed by cell id, values are lists of tuple coordinates.
        tau_range : iterable
            range of time lags to try for calculation.
        interval_range : iterable
            range of interval distances to use for calculation of
            the linregress about point of interest p_n.

        Returns
        ----------
        turn_stats : dict
            keyed by tau, containing dict keyed by interval,
            containing dict keyed by cell, containing proportion of
            time a cell turns right
        theta_stats : dict
            keyed by tau, containing dict keyed by interval,
            containing dict keyed by cell, containing a list of
            three elements [avg_theta, min_theta, max_theta]
        '''
        turn_stats = {}
        theta_stats = {}
        for tau in tau_range:
            turn_stats[tau] = {}
            theta_stats[tau] = {}
            for interval in interval_range:
                turn_stats[tau][interval] = {}
                theta_stats[tau][interval] = {}
                for u in cell_ids:
                    cell_path = cell_ids[u]
                    turns, thetas = self.turn_direction(cell_path, tau = tau, interval = interval)
                    turn_stats[tau][interval][u] = np.mean(turns)
                    # calc theta stats, set as val in triple dict
                    ts = [np.mean(thetas), min(thetas), max(thetas)]
                    theta_stats[tau][interval][u] = ts
        return turn_stats, theta_stats

class MSDFeatures(object):
    '''
    Calculates mean squared displacement related features.

    Attributes
    ----------
    cell_ids : dict
        keyed by cell ids, values are lists of coordinate tuples.
    msd_distributions : dict
        MSDs for each cell across a range of time lags `tau`.
        keyed by cell ids, values are lists of floats.
    log_distributions : dict
        log transform of `msd_distributions`.
        keyed by cell ids, values are lists of floats.
    alphas : dict
        alpha exponent for each cells motion.
        keyed by cell ids, values are floats.

    Notes
    -----
    MSD(tau) = (1 / tau) * sum( (x(t + tau) + x(t))^2 )

    Plotting a distribution of Tau vs MSD will generate a curve.
    The exponential nature of this curve will describe the type
    of motion the particle is experiencing.

    1 < alpha --> impeded diffusion
    Linear curve (alpha = 1) --> diffusion
    1 < alpha < 2 --> super diffusion
    alpha = 2 --> ballistic motion
    '''
    def __init__(self, cell_ids):
        self.msd_distributions, self.log_distributions = self.msd_dist_all_cells(cell_ids, tau = 31)
        self.alphas = self.calc_alphas(self.log_distributions)

    def calc_msd(self, path, tau):
        '''
        Calculates mean squared displacement for a cell path
        for a given time lag tau

        Parameters
        ----------
        path : list
            tuples of sequential XY coordinates.
        tau  : int
            time lag to consider when calculating MSD.

        Returns
        -------
        msd : float
            mean squared displacement of path given time lag tau.
        '''
        distances = []
        t = 0
        while (t + tau) < len(path):
            distsq = ( distance(path[ t + tau ], path[ t ]) )**2
            distances.append(distsq)
            t += 1

        msd = sum(distances)/t

        return msd

    def calc_msd_distribution(self, path, max_tau):
        '''
        Calculates the distribution of MSDs for a range of time
        lag values tau.

        Parameters
        ----------
        path : list
            tuples containing sequential XY coordinates.
        max_tau : int
            maximum time lag `tau` to consider for MSD calculation.

        Returns
        -------
        distribution : list
            MSDs indexed by `tau`.
        log_distribution : list
            log transform of `distribution`.

        Notes
        -----
        If a cell is stationary indefinitely, it effectively has
        a MSD of `0`.
        However, this calls a math domain error, so checks are in
        place that ensure the final `alpha` calculation will be
        0 without raising exceptions.

        Here -- returns both distributions as a string 'flat' if
        MSD calc is 0 for a given range tau
        '''
        tau = 1
        distribution = []
        log_distribution = []
        while tau < max_tau:
            msd = self.calc_msd(path, tau)
            distribution.append(msd)
            if msd > 0:
                log_distribution.append( math.log(msd) )
            else:
                distribution = 'flat'
                log_distribution = 'flat'
                return distribution, log_distribution
            tau += 1

        return distribution, log_distribution

    def msd_dist_all_cells(self, cell_ids, tau = 31):
        '''
        Calculates MSD distributions for all cells in dict cell_ids

        Parameters
        -------------------
        cell_ids : dictionary keyed by cell_id, contains lists of
                   tuples of sequential XY coordinates
        tau      : maximum time lag to calculate for MSD distributions
                   n.b. if path < 30 t's long, uses len(path) as max tau

        Returns
        -------------------
        msd_distributions : dict keyed by cell_id containing list
                            of MSD distributions, indexed by tau
        log_distributions : dict with log transformed lists of
                            MSD distributions

        Notable behavior
        -------------------
        checks for distributions == 'flat' string and sets values for
        that cell as 'flat', which is checked to provide an alpha of 0

        '''

        msd_distributions = {}
        log_distributions = {}
        for u in cell_ids:
            cell = cell_ids[u]
            if len(cell) > 30:
                max_tau = tau
            else:
                max_tau = len(cell)
            dist, log_dist = self.calc_msd_distribution(cell, max_tau)
            if dist == 'flat':
                msd_distributions[u] = 'flat'
                log_distributions[u] = 'flat'
            else:
                msd_distributions[u] = dist
                log_distributions[u] = log_dist

        return msd_distributions, log_distributions

    def calc_alphas(self, log_distributions):
        '''
        Calculate the alpha coefficient value as the slope of a
        log(MSD) v log(tau) plot

        Parameters
        ----------
        log_distributions : dict
            keyed by cell_id containing lists of log transformed MSDs,
            indexed by tau

        Returns
        -------
        alphas : dict
            keyed by cell_id, with values as the alpha coeff
            from log(MSD(tau)) = alpha*log(tau)

        Notes
        -----
        Checkes if log_distributions[cell_id] is == 'flat'
        If so, sets alpha at the appropriate "flat" slope of 0
        '''

        alphas = {}

        for u in log_distributions:
            if log_distributions[u] == 'flat':
                alphas[u] = 0
            else:
                tau = np.arange(1, len(log_distributions[u]) + 1 )
                log_tau = np.log(tau)

                slope, intercept, r_val, p_val, SE = stats.linregress( log_tau, log_distributions[u] )
                alphas[u] = slope

        return alphas

    def _plot_msd_dists(self, output_dir, msd_distributions, log_distributions):
        '''
        Plots tau vs MSD and log(tau) vs log(MSD)
        Saves as a PDF in the motility_statistics output folder
        '''

        tau = range(1,31)
        log_tau = []
        for val in tau:
            log_tau.append( math.log(val) )

        plt.figure(1)
        plt.subplot(121) # call subplots, 1 row, 2 columns, plot 1
        for u in msd_distributions:
            plt.scatter(tau, msd_distributions[u])
            plt.plot(tau, msd_distributions[u])
            plt.xlabel(r'$\tau$') # raw text string, latex to call tau char
            plt.ylabel('Mean Square Displacement')

        plt.subplot(122) # call subplots, 1 row, 2 columns, plot 2
        for u in msd_distributions:
            plt.scatter(log_tau, log_distributions[u])
            plt.plot(log_tau, log_distributions[u])
            plt.xlabel(r'log($\tau$)')
            plt.ylabel('log(Mean Square Displacement)')

        plt.subplots_adjust(wspace = 0.2) # increasing spacing b/w subplots

        plt.savefig(str(output_dir + 'MSD_Plots.pdf'))
        return

class RWFeatures(object):
    '''
    Calculates features relative to a random walk and self-similarity measures.

    Attributes
    ----------
    cell_ids : dict
        keyed by cell ids, values are lists of coordinate tuples.
    gf : hmsims.GeneralFeatures object.
    diff_linearity : dict
        difference in linearity r**2 between observed cell and a random walk.
        keyed by cell ids, values are floats.
    diff_net_dist : dict
        difference in net distance between an observed cell and a random walk.
        keyed by cell ids, values are floats.
    cell_kurtosis : dict
        measured displacement distribution kurtosis.
        keyed by cell ids, values are floats.
    diff_kurtosis : dict
        difference between `cell_kurtosis` and a random walk kurtosis.
        keyed by cell ids, values are floats.
    nongaussalpha : dict
        non-Gaussian parameters alpha_2 of the displacement distribution.
        keyed by cell ids, values are floats.
    disp_var : dict
        variance of the displacement distribution.
        keyed by cell ids, values are floats.
    hurst_RS : dict
        Hurst coefficients of the displacement time series, as estimated by
        Mandelbrot's original Rescaled Range method.
        keyed by cell ids, values are floats.
    autocorr_max_tau : int
        maximum tau for autocorrelation calculation.
    autocorr : dict
        autocorrelation coefficients for displacement time series.
        keyed by cell ids, values are floats.
    '''
    def __init__(self, cell_ids, gf):
        '''
        Calculates features relative to a random walk and self-similarity measures.

        Parameters
        ----------
        cell_ids : dict
            keyed by cell ids, values are lists of coordinate tuples.
        gf : hmsims.GeneralFeatures object.
        '''
        self.cell_ids = cell_ids
        self.gf = gf

        self.diff_linearity, self.diff_net_dist = self.run_comparative_randwalk(cell_ids, gf.linearity, gf.net_distance, gf.avg_speed)
        self.cell_kurtosis, self.diff_kurtosis = self.vartau_kurtosis_comparison(cell_ids, max_tau = 10)

        self.nongaussalpha = self.nongauss_coeff(cell_ids)
        self.disp_var, self.disp_skew = self.displacement_props(cell_ids)
        self.hurst_RS = self.hurst_mandelbrot(cell_ids)
        self.autocorr_max_tau = 11
        self.autocorr, self.qstats, self.pvals = self.autocorr_all(cell_ids, max_tau = self.autocorr_max_tau)

    def random_walk(self, origin, N, speed_mu, speed_sigma = None):
        '''
        Parameters
        ----------
        origin : starting point for the model
        N      : number of steps to model
        speed_mu : mean speed for the model
        speed_sigma : stdev of the normal distribution for step sizes

        Returns
        -------
        model_net_dist : float
        model_linearity : float
        model_path : list

        Notes
        -----
        Net distance of random walks is based on
        tau = time_unit
        rate = sigma = step_size = time_unit * velocity_x
        n = number of steps (can be exchanged for t if =)
        <x(t)^2> = n * sigma^2
        <x(t)^2>^0.5 = sigma * sqrt(n) = step_size * sqrt(n)
        <(x,y)^2> = <r^2> = <x(t)^2> + <y(t)^2> = 2 * n * sigma^2
        <r^2>^0.5 = sqrt(2) * step_size * sqrt(n)

        step_size = rate = c = sqrt(x_step^2 + y_step^2) = sqrt(2*step^2)
        x_step = y_step = rate/sqrt(2)
        (Random Walks in Biology, 1992, Howard C. Berg)
        '''
        model_net_dist = math.sqrt(2) * speed_mu * math.sqrt( N )
        model_path = [ origin ]
        if speed_sigma == None:
            speed_sigma = 0.2*speed_mu
        rate = round( np.random.normal(speed_mu, speed_sigma), 3)
        vectors = [ (0, rate), (0, -rate), (rate, 0), (-rate, 0) ]
        i = 0
        while i < N:
            walk = np.random.random()
            if 0 <= walk < 0.25:
                step = vectors[0]
            elif 0.25 <= walk < 0.50:
                step = vectors[1]
            elif 0.50 <= walk < 0.75:
                step = vectors[2]
            else:
                step = vectors[3]

            new_x = model_path[-1][0] + step[0]
            new_y = model_path[-1][1] + step[1]
            model_path.append( (new_x, new_y) )
            i += 1


        x = []
        y = []
        for coor in model_path:
            x.append( coor[0] )
            y.append( coor[1] )

        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            model_linearity = r_value**2
        except:
            return 'error', 'error', 'error'

        return model_net_dist, model_linearity, model_path

    def compare_to_randwalk(self, cell, u, cell_linearity, cell_net_dist, cell_avg_speed):
        '''
        Compares an observed path to a simulated random walk with the
        same displacement mean
        '''
        origin = cell[0]
        N = len(cell)
        rate = cell_avg_speed
        linearitys = []

        i = 0
        while i < 200:
            model_net_dist, model_linearity, model_path = self.random_walk(origin, N, rate)
            if model_net_dist == 'error':
                pass
            else:
                linearitys.append( model_linearity )
            i += 1

        avg_linearity = np.mean(linearitys)

        diff_linearity = cell_linearity - avg_linearity
        diff_net_dist = cell_net_dist - model_net_dist

        return diff_linearity, diff_net_dist

    def run_comparative_randwalk(self, cell_ids, linearity, net_distance, avg_speed):
        '''Runs comparisons to random walks for all cells in `cell_ids`'''
        diff_linearity = {}
        diff_net_dist = {}
        for u in cell_ids:
            cell = cell_ids[u]
            lin, net = self.compare_to_randwalk(cell, u, linearity[u], net_distance[u], avg_speed[u])
            diff_linearity[u] = lin
            diff_net_dist[u] = net

        return diff_linearity, diff_net_dist

    def make_displacement_dist(self, cell_ids, tau):
        '''Generates displacement distributions for all cells in `cell_ids`'''

        displacement_dist = {}
        for u in cell_ids:
            cell = cell_ids[u]
            displacement_dist[u] = []
            i = 1
            while (i + tau) < len(cell):
                d = distance( cell[i + tau], cell[i] )
                displacement_dist[u].append( d )
                i = i + tau

        return displacement_dist

    def kurtosis_comparison(self, cell_ids, displacements):
        '''Compares cell kurtosis to random walk kurtosis'''
        cell_kurtosis = {}
        diff_kurtosis = {}
        rayleigh_kurtosis = 3.245089
        for u in cell_ids:
            dist = displacements[u]
            cell_kurtosis[u] = stats.kurtosis(dist)
            diff_kurtosis[u] = cell_kurtosis[u]

        return cell_kurtosis, diff_kurtosis

    def vartau_kurtosis_comparison(self, cell_ids, max_tau = 10 ):
        '''
        Calculate `cell_kurtosis` and `diff_kurtosis` for a range of possible
        time intervals `tau`.
        '''
        # set up dicts keyed by tau, containing dicts keyed by cell_id
        all_cell_kurtosis = {}
        all_diff_kurtosis = {}

        tau = 1
        while tau <= max_tau:
            displacement_dist = self.make_displacement_dist(cell_ids, tau)
            all_cell_kurtosis[tau], all_diff_kurtosis[tau] = self.kurtosis_comparison(cell_ids, displacements = displacement_dist)
            tau += 1

        return all_cell_kurtosis, all_diff_kurtosis

    def moving_kurt_comparison(self, cell_ids, displacements, speed):
        '''Compare kurtosis only while cell is moving'''
        cell_kurtosis = {}
        diff_kurtosis = {}
        rayleigh_kurtosis = 3.245089
        for u in cell_ids:
            dist = np.array(displacements[u])
            dist = dist[dist > speed]
            if len(dist) > 0:
                cell_kurtosis[u] = stats.kurtosis(dist)
                diff_kurtosis[u] = cell_kurtosis[u] - rayleigh_kurtosis
            else:
                cell_kurtosis[u] = 0
                diff_kurtosis[u] = 0

        return cell_kurtosis, diff_kurtosis

    def moving_kurt(self, cell_ids, max_tau = 10, moving_range = range(1,11)):
        '''
        Parameters
        ----------

        cell_ids : dict
            keyed by cell_id, values are lists of corrdinate tuples.
        max_tau : int
            maximum time lag to consider for kurtosis calculations.
        moving_range : iterable
            range of speeds to consider as a threshold for movement.

        Returns
        -------
        cell_moving_kurt : triple dict raw kurtosis of all cells in cell_ids
                           keyed by speed threshold, tau, and cell_id
                           cell_moving_kurt = {
                            speed1 : {
                                tau1 :
                                    {cell_id1 : kurt
                                    ...}
                                }
                                ...
                            } ...
                           }
        diff_moving_kurt : kurtosis of all cells, normalized by Rayleigh kurt
        structured as above
        '''

        cell_moving_kurt = {}
        diff_moving_kurt = {}

        for speed in moving_range:
            cell_moving_kurt[speed] = {}
            diff_moving_kurt[speed] = {}
            tau = 1
            while tau <= max_tau:
                displacement_dist = self.make_displacement_dist(cell_ids, tau)
                cell_moving_kurt[speed][tau], diff_moving_kurt[speed][tau] = self.moving_kurt_comparison(cell_ids, displacements= displacement_dist, speed=speed)

        return cell_moving_kurt, diff_moving_kurt

    def displacement_props(self, cell_ids):
        '''
        Calculates variance and skewness of the displacement distribution
        for each cell in cell_ids

        Parameters
        ----------
        cell_ids : dict of lists containing tuples of sequential XY coordinates

        Returns
        -------
        var : dict keyed by cell_id with variance of displacement distribution
        skew : dict keyed by cell_ids with skew of displacement distribution
        '''
        var = {}
        skew = {}

        allX = self.make_displacement_dist(cell_ids=cell_ids, tau = 1)
        for u in allX:
            X = np.asarray( allX[u] )
            var[u] = np.var(X)
            skew[u] = stats.skew(X)

        return var, skew

    def calc_ngaussalpha(self, X):
        '''
        Calculates the non-Gaussian parameter alpha_2 of a given displacement
        distribution, X

        Parameters
        ----------
        X : array-like
            distribution of displacements.

        Returns
        -------
        alpha_2 : float
            non-Gaussian coefficient, floating point values.

        Notes
        -----
        alpha_2 = <dx^4> / 3*<dx^2>^2 - 1

        effectively, a ratio of coeff* (kurtosis / variance)

        For a Gaussian distribution, alpha_2 = 0
        For non-Gaussian distributions, alpha_2 increases with the length
        of the tails

        Levy-like motion would be expected to have alpha_2 > 0, while diffusive
        motion would be expected to have alpha_2 == 0

        Rayleigh alpha_2 ~= -0.33
        '''

        X = np.asarray(X)
        alpha_2 = np.mean( X**4 ) / (3*np.mean(X**2)**2) - 1

        return alpha_2

    def nongauss_coeff(self, cell_ids):
        '''
        Calculates non-Gaussian coefficient alpha_2 for all cells in cell_ids.

        Parameters
        ----------
        cell_ids : dict
            lists containing tuples of sequential XY coordinates.

        Returns
        -------
        nongauss_coeff : dict
            keyed by cell_id, nonGauss coeff values.
        '''

        ngaussalpha = {}
        allX = self.make_displacement_dist(cell_ids=cell_ids, tau=1)
        for u in allX:
            X = allX[u]
            ngaussalpha[u] = self.calc_ngaussalpha(X)

        return ngaussalpha


    def largest_pow2(self, num):
        '''Finds argmax_n 2**n < `num`'''
        for i in range(0,10):
            if int(num / 2**i) > 1:
                continue
            else:
                return i-1

    def rescaled_range(self, X, n):
        '''
        Finds rescaled range <R(n)/S(n)> for all sub-series of size `n` in `X`
        '''
        N = len(X)
        if n > N:
            return None
        # Create subseries of size n
        num_subseries = int(N/n)
        Xs = np.zeros((num_subseries, n))
        for i in range(0, num_subseries):
            Xs[i,] = X[ int(i*n) : int(n+(i*n)) ]

        # Calculate mean rescaled range R/S
        # for subseries size n
        RS = []
        for subX in Xs:

            m = np.mean(subX)
            Y = subX - m
            Z = np.cumsum(Y)
            R = max(Z) - min(Z)
            S = np.std(subX)
            if S <= 0:
                print("S = ", S)
                continue
            RS.append( R/S )
        RSavg = np.mean(RS)

        return RSavg

    def hurst_mandelbrot(self, cell_ids):
        '''
        Calculates the Hurst coefficient `H` for each cell in `cell_ids`.

        Notes
        -----
        for E[R(n)/S(n)] = Cn**H as n --> inf
        H : 0.5 - 1 ; long-term positive autocorrelation
        H : 0.5 ; fractal Brownian motion
        H : 0-0.5 ; long-term negative autocorrelation

        N.B. SEQUENCES MUST BE >= 18 units long, otherwise
        linear regression for log(R/S) vs log(n) will have
        < 3 points and cannot be performed
        '''

        hurst_RS = {}
        allX = self.make_displacement_dist(cell_ids, tau = 1)
        for u in allX:
            X = allX[u]
            RSl = []
            ns = []
            for i in range(0,self.largest_pow2(len(X))):
                ns.append(int(len(X)/2**i))
            for n in ns:
                RSl.append( self.rescaled_range(X, n) )
            m, b, r, pval, sderr = stats.linregress(np.log(ns), np.log(RSl))
            hurst_RS[u] = m
        return hurst_RS

    # using Dominik Krzeminski's (@dokato) dfa library
    # see
    # https://github.com/dokato/
    # https://github.com/dokato/dfa/blob/master/dfa.py
    def dfa_all(self, cell_ids):
        '''
        Performs detrended fluctuation analysis on cell displacement series.

        Parameters
        ----------
        cell_ids : dict of lists keyed by cell_id
        ea. list represents a cell. lists contain sequential tuples
        containing XY coordinates of a cell at a given timepoint

        Returns
        -------
        dfa_alpha : dict keyed by cell_id, values are the alpha coefficient
        calculated by detrended fluctuation analysis

        References
        ----------
        Ping et. al., Mosaic organization of DNA nucleotides, 1994, Phys Rev E
        '''
        dfa_alpha = {}
        allX = self.make_displacement_dist(cell_ids, tau = 1)
        for u in allX:
            X = allX[u]
            scales, fluct, coeff = dfa(X, scale_lim = [2,7], scale_dens = 0.25)
            dfa_alpha[u] = coeff
        return dfa_alpha

    def autocorr_all(self, cell_ids, max_tau = 10):
        '''
        Estimates the autocorrelation coefficient for each series of cell
        displacements over a range of time lags.

        Parameters
        ----------
        cell_ids : dict of lists keyed by cell_id
        ea. list represents a cell. lists contain sequential tuples
        containing XY coordinates of a cell at a given timepoint

        Returns
        -------
        autocorr : dict of lists, containing autocorrelation coeffs for
        sequential time lags
        qstats : dict of lists containing Q-Statistics (Ljung-Box)
        pvals : dict of lists containing p-vals, as calculated from Q-Statistics

        Notes
        -----
        Estimation method:
        https://en.wikipedia.org/wiki/Autocorrelation#Estimation

        R(tau) = 1/(n-tau)*sigma**2 [sum(X_t - mu)*(X_t+tau - mu)] | t = [1,n-tau]

        X as a time series, mu as the mean of X, sigma**2 as variance of X
        tau as a given time lag (sometimes referred to as k in literature)

        Implementation uses statsmodels.tsa.stattools.acf()

        n.b. truncated to taus [1,10], to expand to more time lags, simply
        alter the indexing being loaded into the return dicts
        '''
        autocorr = {}
        qstats = {}
        pvals = {}

        allX = self.make_displacement_dist(cell_ids, tau = 1)
        for u in allX:
            X = allX[u]
            # Perform autocorrelation, lags [0,29] tested
            # Perform Ljung-Box Q-statistic calculation to determine if
            # autocorrelations detected are significant or random
            ac, q, p = acf(X, adjusted = True, nlags = (max_tau+1), qstat=True)
            autocorr[u] = ac[1:max_tau]
            qstats[u] = q[1:max_tau]
            pvals[u] = p[1:max_tau]

        return autocorr, qstats, pvals

    def partial_acorr_all(self, cell_ids, max_tau = 10):
        '''
        Estimates the partial autocorrelation coefficient for each series of cell
        displacements over a range of time lags.

        Parameters
        ----------
        cell_ids : dict
            keyed by cell_id, values are lists of coordinate tuples.
        max_tau : int
            maximum tau to estimate partial autocorrelation.

        Returns
        -------
        partial_acorr : dict
            keyed by cell id, values are lists, containing autocorrelation coeffs
            for sequential time lags

        Notes
        -----
        n.b. truncated to taus [1,10], to expand to more time lags, simply
        alter the indexing being loaded into the return dicts.
        Implementation uses statsmodels.tsa.stattools.pacf()
        '''
        partial_acorr = {}

        allX = self.make_displacement_dist(cell_ids, tau = 1)
        for u in allX:
            X = allX[u]
            n = len(X)
            pac = pacf(X, nlags = max_tau)
            partial_acorr[u] = pac[1:max_tau]

        return partial_acorr
