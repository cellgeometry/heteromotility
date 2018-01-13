from __future__ import division
import random
from .hmstats import GeneralFeatures, RWFeatures, MSDFeatures
from .hmtrack import CellPaths
from .hmtools import dictofdict2array, tripledict2array, merge_flat_lists
import numpy as np
from functools import partial
from multiprocessing import Pool


'''
#---------------------------
# MODULE CONTENTS
#---------------------------
Classes to generate simulated data sets for analysis, including
1. Random walk simulations
2. Biased random walk simulations
3. Levy and power flight simulations

#---------------------------
# To Do
#---------------------------

1. Parallelize simulations with multiprocess.Pool.map() and functools.partial()

'''


class RWSims(object):
    def __init__(self, obj = 100, length = 100, bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf), speed_range=range(5,10), bias_prop=0.5):
        self.rw_paths = self.rwsim(obj, length, bound_x=bound_x, bound_y=bound_y, speed_range=speed_range)
        self.brw_paths = self.biased_rwsim(obj, length, bound_x=bound_x, bound_y=bound_y, speed_range=speed_range, bias_prop=bias_prop)

    def random_walk(self, N, speed_mu = 7.0, speed_sigma = None, origin = (0,0), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf)):
        import numpy as np
        '''
        Models random walk given an origin, length N, step size : rate
        abbreviated from implementation in hmstats.RWFeatures

        Parameters
        -------------
        N : length of walk to model
        speed_mu : float, optional.
            mean of Gaussian distribution speeds are sampled from
        speed_sigma : float, optional.
            standard deviation of Gaussian distribution speeds are sampled from
        origin : tuple.
            origin point, default = (0,0)
        bound_x : tuple, optional.
            (min, max) allowed values for X coordinates.
            Useful to confine model to a specified box.
        bound_y : tuple, optional.
            (min, max) allowed values for Y coordinates.
            Useful to confine model to a specified box.

        Returns
        -------------
        model_path : list containing sequential tuples of XY coordinates

        Notes
        -------------
        step_size = rate = c = sqrt(x_step^2 + y_step^2) = sqrt(2*step^2)
        x_step = y_step = rate/sqrt(2)
        (Random Walks in Biology, 1992, Howard C. Berg)
        '''

        model_path = [ origin ]
        if speed_sigma == None:
            speed_sigma = 0.2*speed_mu
        i = 0
        while i < N:
            rate = round(np.random.normal(speed_mu, speed_sigma), 3)

            if i == 0:
                direction = np.random.random(1) * 2 * np.pi

            turn = np.random.random() * 2 * np.pi
            direction = (direction + turn) % (2 * np.pi)

            dx = np.cos(direction)*rate
            dy = np.sin(direction)*rate

            step = (float(dx), float(dy))

            # vectors = [ (0, rate), (0, -rate), (rate, 0), (-rate, 0) ]
            # walk = np.random.random()
            # if 0 <= walk < 0.25:
            #     step = vectors[0]
            # elif 0.25 <= walk < 0.50:
            #     step = vectors[1]
            # elif 0.50 <= walk < 0.75:
            #     step = vectors[2]
            # else:
            #     step = vectors[3]

            try_x = model_path[-1][0] + step[0]
            try_y = model_path[-1][1] + step[1]
            if bound_x[0] < try_x < bound_x[1] and bound_y[0] < try_y < bound_y[1]:
                model_path.append( (try_x, try_y) )
                i += 1
            else:
                pass


        return model_path

    def rwsim(self, obj, length, speed_range = range(5, 16), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf)):
        '''
        Generates simulated object paths using an unbiased random walk model

        Parameters
        ------------
        obj : integer specifying the number of object paths to be simulated
        length : integer specifying the length of each path to be simulated
        speed_range : tuple containing integers of step sizes to use when
        building models
        bound_x : tuple, optional.
            (min, max) allowed values for X coordinates.
            Useful to confine model to a specified box.
        bound_y : tuple, optional.
            (min, max) allowed values for Y coordinates.
            Useful to confine model to a specified box.

        Returns
        ------------
        rw_paths : dict keyed by speed, containing dicts keyed by obj id,
        containing lists of tuples, each tuple
        specifying the XY position of a simulated object at a given timepoint
            rw_paths =
            {
            speed_i : {0: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
                        1: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
                        ...},
            speed_i+1 : ...,
            ...
            }
        '''

        rw_paths = {}

        for speed in speed_range:
            rw_paths[speed] = {}

            '''
            rwfun = partial(self.random_walk, N = length, speed_mu = speed, origin = (0,0), bound_x=bound_x, bound_y=bound_y)
            def mp_wrapper(_):
                return rwfun()

            p = Pool()
            model_paths = p.map(mp_wrapper, range(obj))

            for i in range(obj):
                rw_paths[speed][i] = model_paths[i]
            '''

            i = 0
            while i < obj:
                model_path = self.random_walk(N = length, speed_mu = speed, origin = (0,0), bound_x=bound_x, bound_y=bound_y)
                rw_paths[speed][i] = model_path
                i += 1

        return rw_paths

    def biased_random_walk(self, N, speed_mu = 7.0, speed_sigma = None, bias_prop = 0.5 , origin = (0,0), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf)):
        import numpy as np
        '''
        * Models a set of random walks biased to move in one direction
        * Randomly selects in which direction to bias movement
        * Object will move in the direction of bias a proportion of the time
        specified in bias_prop

        Parameters
        ------------
        N : length of walk to be modeled
        rate : step size to be modeled
        bias_prop : proportion of the time an object moves in the direction of bias
        origin : origin point for the start of the walk to be modeled
        bound_x : tuple, optional.
            (min, max) allowed values for X coordinates.
            Useful to confine model to a specified box.
        bound_y : tuple, optional.
            (min, max) allowed values for Y coordinates.
            Useful to confine model to a specified box.
        Returns
        ------------
        model_path : list of tuples, each containing an XY point

        '''
        model_path = [origin]
        if speed_sigma == None:
            speed_sigma = 0.2*speed_mu
        bias_choice = random.randrange(start = 0, stop = 4, step = 1)
        bias_direction = np.random.random(1) * 2 * np.pi

        i = 0
        while i < N:

            if i == 0:
                direction = np.random.random(1) * 2 * np.pi

            if np.random.random() > bias_prop:
                turn = np.random.random() * 2 * np.pi
                direction = (direction + turn) % (2 * np.pi)
            else:
                direction = bias_direction

            rate = round(np.random.normal(speed_mu, speed_sigma), 3)

            dx = np.cos(direction)*rate
            dy = np.sin(direction)*rate

            step = (float(dx), float(dy))

            # walk = np.random.random()
            # vectors = [ (0, rate), (0, -rate), (rate, 0), (-rate, 0) ]
            # bias = vectors[bias_choice]
            # vectors.pop(bias_choice)
            # unbiased_prop = 1-bias_prop
            # break1 = bias_prop + (unbiased_prop/3)
            # break2 = bias_prop + 2*(unbiased_prop/3)
            # if walk <= bias_prop:
            #     step = bias
            # elif bias_prop < walk <= break1:
            #     step = vectors[0]
            # elif break1 < walk <= break2:
            #     step = vectors[1]
            # else:
            #     step = vectors[2]

            try_x = model_path[-1][0] + step[0]
            try_y = model_path[-1][1] + step[1]
            if bound_x[0] < try_x < bound_x[1] and bound_y[0] < try_y < bound_y[1]:
                model_path.append( (try_x, try_y) )
                i += 1
            else:
                pass


        return model_path

    def biased_rwsim(self, obj, length, speed_range = range(5,16), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf), bias_prop=0.5):
        '''
        Generates simulated object paths using a biased random walk model

        Parameters
        ------------
        obj : integer specifying the number of object paths to be simulated
        length : integer specifying the length of each path to be simulated

        Returns
        ------------
        brw_paths : dict keyed by obj id, containing lists of tuples, each tuple
        specifying the XY position of a simulated object at a given timepoint
            rw_paths = {
            0: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            1: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            ...
            }
        '''

        brw_paths = {}

        for speed in speed_range:
            brw_paths[speed] = {}
            i = 0
            while i < obj:
                model_path = self.biased_random_walk(N = length, speed_mu = speed, bias_prop = bias_prop, origin = (0,0), bound_x=bound_x, bound_y=bound_y)
                brw_paths[speed][i] = model_path
                i += 1

        return brw_paths

class LevySims(object):
    def __init__(self, obj = 100, length = 100,
                bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf),
                p_scale = range(25,30), alpha = 2.000, levy=False, levy_scale = range(1,5)):
        '''
        Simulated optimal power and true 'Levy' flier walker models.

        Parameters
        ----------
        obj : integer. number of objects to simulate.
        length : integer. length of simulated tracks.
        bound_x : tuple, (min, max). bounding values on the x axis.
        bound_y : tuple, (min, max). bounding values on the y axis.
        p_scale : iterable. scaling values for power fliers. 25 = 5 px step mean.
        alpha : float. alpha exponent for power distribution.
        levy : boolean. simulate Levy fliers.
        levy_scale : iterable. scaling values for Levy fliers.
        '''
        self.power_paths = self.power_sim(obj, length, bound_x=bound_x, bound_y=bound_y, scale_range = p_scale, alpha=alpha)
        if levy:
            self.levy_paths = self.levy_sim(obj, length, bound_x=bound_x, bound_y=bound_y, scale_range = levy_scale)

    def power_flier( self, N, x_min = 0.1, alpha = 2.000, scale_coeff=25, origin = (0,0), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf) ):
        '''
        Defines a levy flight path beginning at the origin, taking N steps,
        with step sizes pulled from a power law distribution

        Directions are chosen randomly, as in a random walk
        Step sizes however are chosen according to the distribution:

            P(l_j) ~ l_j**(-a)
            where a is a parameter in the range 1 < u <= 3

        Defaults to a = 2 as the optimal parameter for foraging behavior
        see : Viswanathan 1999, Nature

        Obtaining Levy distributed random variables
        -------------
        To transform uniformly distributed random variable
        into another distribution, use inverse cum. dist. func. (CDF)

        if F is CDF corresponding to probability density of variable f and u
        is a uniformly distributed random variable on [0,1]

        x = F^-1(u)
        for pure power law distribution P(l) = l**-a
        F(x) = 1 - ( x / x_min )
        where x_min is a lower bound on the random variable

        F^-1(u) = x_min(1 - u)^(-1/a)

        See: (Devroye 1986, Non-uniform random variable generation)

        or

        X = F^-1(U) = c / ( phi^-1(1-U/2) )**2 + mu

        c : scaling parameter
        Phi(x) : CDF of Gaussian distribution
        mu : location

        for a pure Levy distribution

        see http://www.math.uah.edu/stat/special/Levy.html

        Parameters
        -------------
        N : number of steps in the flight
        step_min : integer value, minimum step size
        origin : initial site of walk
        u = parameter for probability distribution of step sizes

        Returns
        -------------
        model_path : a list containing tuples of XY coordinates
        '''

        model_path = [ origin ]
        i = 0
        while i < N:

            rate = ( x_min*( 1-random.random() )**(-1/alpha) ) * scale_coeff

            if i == 0:
                direction = np.random.random(1) * 2 * np.pi

            turn = np.random.random() * 2 * np.pi
            direction = (direction + turn) % (2 * np.pi)

            dx = np.cos(direction)*rate
            dy = np.sin(direction)*rate

            step = (float(dx), float(dy))

            try_x = model_path[-1][0] + step[0]
            try_y = model_path[-1][1] + step[1]
            if bound_x[0] < try_x < bound_x[1] and bound_y[0] < try_y < bound_y[1]:
                model_path.append( (try_x, try_y) )
                i += 1
            else:
                pass

        return model_path

    def power_sim(self, obj, length, scale_range = range(25,35), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf), alpha = 2.000):
        '''
        Generates simulated object paths using a levy flight model with
        a power law distribution of displacements

        Parameters
        ------------
        obj : integer specifying the number of object paths to be simulated
        length : integer specifying the length of each path to be simulated

        Returns
        ------------
        levy_paths : dict keyed by obj id, containing lists of tuples, each tuple
        specifying the XY position of a simulated object at a given timepoint
            levy_paths = {
            0: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            1: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            ...
            }

        Distribution with scale factor = 25 has mean = 5
        '''

        power_paths = {}

        for scale in scale_range:
            power_paths[scale] = {}
            i = 0
            while i < obj:
                model_path = self.power_flier(N = length, scale_coeff=scale, bound_x=bound_x, bound_y=bound_y, alpha = alpha)
                power_paths[scale][i] = model_path
                i += 1

        return power_paths

    def levy_flier(self, N, scale_coeff = 1, origin = (0,0), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf) ):
        '''
        Defines a levy flight path beginning at the origin, taking N steps,
        with step sizes pulled from a true Levy distribution

        scale = 2, median disp = 4.4
        scale = ~2.3 has a median displacement of 5 units

        Directions are chosen randomly, as in a random walk
        Step sizes however are chosen according to the distribution:

            P(l_j) ~ l_j**(-a)
            where a is a parameter in the range 1 < u <= 3

        Defaults to a = 2 as the optimal parameter for foraging behavior
        see : Viswanathan 1999, Nature

        Obtaining Levy distributed random variables
        -------------
        To transform uniformly distributed random variable
        into another distribution, use inverse cum. dist. func. (CDF)

        if F is CDF corresponding to probability density of variable f and u
        is a uniformly distributed random variable on [0,1]

        x = F^-1(u)
        for pure power law distribution P(l) = l**-a
        F(x) = 1 - ( x / x_min )
        where x_min is a lower bound on the random variable

        F^-1(u) = x_min(1 - u)^(-1/a)

        See: (Devroye 1986, Non-uniform random variable generation)

        or

        X = F^-1(U) = c / ( phi^-1(1-U/2) )**2 + mu

        c : scaling parameter
        Phi(x) : CDF of Gaussian distribution
        mu : location

        for a pure Levy distribution

        see http://www.math.uah.edu/stat/special/Levy.html

        Parameters
        -------------
        N : number of steps in the flight
        step_min : integer value, minimum step size
        origin : initial site of walk
        u = parameter for probability distribution of step sizes

        Returns
        -------------
        model_path : a list containing tuples of XY coordinates
        '''
        from scipy.stats import levy

        model_path = [ origin ]
        i = 0
        while i < N:

            rate = levy.rvs(scale=scale_coeff)

            if i == 0:
                direction = np.random.random(1) * 2 * np.pi

            turn = np.random.random() * 2 * np.pi
            direction = (direction + turn) % (2 * np.pi)

            dx = np.cos(direction)*rate
            dy = np.sin(direction)*rate

            step = (float(dx), float(dy))

            # vectors = [ (0, rate), (0, -rate), (rate, 0), (-rate, 0) ]
            # walk = random.random()
            # if 0 <= walk < 0.25:
            #     step = vectors[0]
            # elif 0.25 <= walk < 0.50:
            #     step = vectors[1]
            # elif 0.50 <= walk < 0.75:
            #     step = vectors[2]
            # else:
            #     step = vectors[3]

            try_x = model_path[-1][0] + step[0]
            try_y = model_path[-1][1] + step[1]
            if bound_x[0] < try_x < bound_x[1] and bound_y[0] < try_y < bound_y[1]:
                model_path.append( (try_x, try_y) )
                i += 1
            else:
                pass

        return model_path

    def levy_sim(self, obj, length, scale_range = range(1,10), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf) ):
        '''
        Generates simulated object paths using a levy flight model with a true
        Levy distribution of displacements

        Parameters
        ------------
        obj : integer specifying the number of object paths to be simulated
        length : integer specifying the length of each path to be simulated

        Returns
        ------------
        levy_paths : dict keyed by obj id, containing lists of tuples, each tuple
        specifying the XY position of a simulated object at a given timepoint
            levy_paths = {
            0: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            1: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            ...
            }
        '''

        levy_paths = {}

        for scale in scale_range:
            levy_paths[scale] = {}
            i = 0
            while i < obj:
                model_path = self.levy_flier(N = length, scale_coeff=scale, bound_x=bound_x, bound_y=bound_y)
                levy_paths[scale][i] = model_path
                i += 1

        return levy_paths


class fBmsims(object):
    def __init__(self, obj = 100, length = 100, H=0.8, bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf), scale_range=range(1,5) ):
        self.fbm_paths = self.fbm_sim(obj, length, H, scale_range = scale_range, bound_x=bound_x, bound_y=bound_y)


    def fbm_mat( self, gridN = 100, T = 100, H = 0.99 ):
        '''
        Generates a Cholesky decomposition of a covariance matrix

        Gamma = (R(i,j), i, j) for i,j = 0...n

        where n is a grid size parameter

        Parameters
        -------------
        gridN : size of grid to generate for the covariance matrix
        T : length of displacement series to generate
        H : hurst parameter [0,1]

        Returns
        -------------
        X : series of displacements following fBm for the given Hurst index
        '''

        covr = np.arange(0,30)
        tcol = np.matrix(np.linspace(0,gridN,T)[1:])
        trow = tcol.T
        Gamma = 0.5*(np.power(trow, 2*H) + np.power(tcol, 2*H) - np.power(np.abs((trow - tcol)), 2*H) )
        Gamma_c = np.linalg.cholesky(Gamma)
        return Gamma_c

    def fbm_series(self, gridN = 100, T = 100, H = 0.99):
        '''
        Generates a series of displacements simulating a fractal
        Brownian process with a given Hurst parameter of length T using
        the Cholesky decomposition method, as outlined in

        T. Dieker, Simulation of fractional Brownian motion, 2004

        Briefly:
        1. Generate a matrix Gamma = [ R(t_i,t_j), for i,j = 0...n ]
        2. take the cholesky decomposition of Gamma
        3. Find Sigma = cholesky_Gamme**0.5
        4. Generate a vector v or length n with Gaussian distributed values
        5. vector u = Sigma * v represents a fBm series
        '''

        sigma = self.fbm_mat(gridN = gridN, T = T, H = H)
        v = np.random.normal(size = sigma.shape[1])
        v = v.reshape(v.size, 1)
        u = sigma*v

        return u

    def fbm_walker(self, N, scale=1, H=0.99, origin=(0,0), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf)):

        model_path = [ origin ]
        i = 0
        u = self.fbm_series(gridN = N+1, T = N+1, H = H)
        while i < N:
            rate = float(u[i])*scale

            if i == 0:
                direction = np.random.random(1) * 2 * np.pi

            turn = np.random.random() * 2 * np.pi
            direction = (direction + turn) % (2 * np.pi)

            dx = np.cos(direction)*rate
            dy = np.sin(direction)*rate

            step = (float(dx), float(dy))

            try_x = model_path[-1][0] + step[0]
            try_y = model_path[-1][1] + step[1]
            if bound_x[0] < try_x < bound_x[1] and bound_y[0] < try_y < bound_y[1]:
                model_path.append( (try_x, try_y) )
                i += 1
            else:
                pass

        return model_path

    def fbm_sim(self, obj, length, H, scale_range = np.linspace(0.1, 0.4, 10), bound_x=(-np.Inf, np.Inf), bound_y=(-np.Inf, np.Inf)):
        '''
        Generates simulated object paths using a fractal Brownian motion model

        Parameters
        ------------
        obj : integer specifying the number of object paths to be simulated
        length : integer specifying the length of each path to be simulated
        H : Hurst coefficient for generation of fBm displacement series
        scale_range : range of scales to use

        Returns
        ------------
        fbm_paths : dict keyed by scale, containing dict keyed by obj id,
        containing lists of tuples, each tuple
        specifying the XY position of a simulated object at a given timepoint
            fbm_paths = {
            1: {
            0: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            1: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            ...},
            2: {
            0: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            1: [ (x1,y1), (x2,y2)...(x_n, y_n) ]
            ...}
            }

        scale = 0.277 is equivalent to mean displacement = 5 units
        '''
        fbm_paths = {}

        for scale in scale_range:
            fbm_paths[scale] = {}
            i = 0
            while i < obj:
                model_path = self.fbm_walker(N = length, scale = scale, H = H, bound_x = bound_x, bound_y = bound_y)
                fbm_paths[scale][i] = model_path
                i += 1

        return fbm_paths

class MultiSim(object):
    def __init__(self, obj = 1000, tau = 20):
        self.pwr2fbm_paths = self.power_fbm(obj = obj, tau = tau)
        self.pwr2urw_paths = self.power_urw(obj = obj, tau = tau)
        self.fbm2urw_paths = self.fbm_rw(obj = obj, tau = tau)

    def join_paths(self, ids1, ids2):
        ids_joined = {}
        for u in ids1:
            ids_joined[u] = ids1[u] + ids2[u]
        return ids_joined

    def power_fbm(self, obj = 1000, tau = 20):
        lfs = LevySims(obj = obj, length = tau, p_scale=range(59,61))
        fbms = fBmsims(obj = 1, length = 1)
        pwr_paths = lfs.power_paths[59]

        comb_paths = {}
        for u in pwr_paths:
            fbm_paths = fbms.fbm_walker(N = tau, scale=1, H=0.99, origin= pwr_paths[u][-1])
            comb_paths[u] = pwr_paths[u] + fbm_paths

        return comb_paths

    def power_urw(self, obj = 1000, tau = 20):
        lfs = LevySims(obj = obj, length = tau, p_scale=range(59,61))
        rws = RWSims(obj = 1, length = 1)
        pwr_paths = lfs.power_paths[59]

        comb_paths = {}
        for u in pwr_paths:
            urw_paths = rws.random_walk(N = tau, speed_mu = 5, origin = pwr_paths[u][-1])
            comb_paths[u] = pwr_paths[u] + urw_paths

        return comb_paths

    def fbm_rw(self, obj = 1000, tau = 20):
        fbms = fBmsims(obj = obj, length = tau)
        rws = RWSims(obj = 1, length = 1)
        fbm_paths = fbms.fbm_paths[4]

        comb_paths = {}
        for u in fbm_paths:
            urw_paths = rws.random_walk(N = tau, speed_mu = 5, origin = fbm_paths[u][-1])
            comb_paths[u] = fbm_paths[u] + urw_paths

        return comb_paths




#---------------------------
# Export Simulated Data
#---------------------------
'''
rws = RWSims(obj = 1000, length = 100)
lfs = LevySims(obj = 1000, length = 100)

for sim in [rws, lfs]:
    cp = CellPaths(cell_ids = rws)
    gf = GeneralFeatures(cell_ids = cp.cell_ids)
    msdf = MSDFeatures(cell_ids = cp.cell_ids)
    rwf = RWFeatures(cell_ids = cp.cell_ids, gf = gf)

    ind_outputs = single_outputs_list(cp.cell_ids, gf = gf, rwf = rwf, msdf = msdf, 'sim_output')
    diff_kurtosis_array = dictofdict2array(rwf.diff_kurtosis)
    avg_moving_speed_array = dictofdict2array( gf.avg_moving_speed )
    time_moving_array = dictofdict2array( gf.time_moving )
    turn_list = tripledict2array(gf.turn_stats)
    theta_list = tripledict2array(gf.theta_stats)
    merged_list = merge_flat_lists([ind_outputs, diff_kurtosis_array, avg_moving_speed_array, time_moving_array, turn_list, theta_list])
'''
