# -*- coding: utf-8 -*-
'''

    The MIT License (MIT)

    Copyright (c) 2014 Kyle Huston

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.

    =============================================================================

    Written by Kyle Huston at 2:39 AM Tuesday, July 8, 2014

    Algorithm based on
        Kaufmann, Bernhard. “Fitting a Sum of Exponentials to Numerical Data.”
        arXiv:physics/0305019, May 6, 2003.
        http://arxiv.org/abs/physics/0305019.
'''

import itertools
import numpy as np
import scipy.optimize as opt
from scipy.misc import factorial

def Cmax(N):
    return int(max([factorial(N)/(factorial(i)*factorial(N-i)) for i in range(N+1)]))

def ComputeAlphaBeta(N,a,b):
    alpha = np.zeros((N+1,Cmax(N)+1))
    beta = np.zeros((N+1,Cmax(N)+1))
    sum_a = np.sum(a)
    beta[0][1] = 1
    alpha[0][1] = sum_a
    for i in np.arange(1,N+1):
        j = 1
        for combo in itertools.combinations(np.arange(N),i):
            beta[i][j] = 1
            alpha[i][j] = sum_a
            for b_index in combo:
                beta[i][j] = beta[i][j] * b[b_index]
                alpha[i][j] = alpha[i][j] - a[b_index+1]
            j = j + 1
    return alpha,beta

def ComputeI0FromParams(x,a,b):
    a = np.array(a)
    b = np.array(b)
    return np.sum(a[1:]*np.exp(-np.outer(b,x)).T,axis=1)

def ComputeIkFromParams(N,x,a,b):
    I0 = ComputeI0FromParams(x,a,b)
    I = []
    I.append(I0)
    if N == 0:
        return I
    else:
        for k in range(N):
            Ik = []
            for t in range(len(x)):
                Ik.append(np.trapz(I[-1][:t+1],x[:t+1]))
            I.append(Ik)
        return np.array(I)

def ComputeI0FromY(x,y):
    return y

def ComputeIkFromY(N,x,y):
    I0 = ComputeI0FromY(x,y)
    I = []
    I.append(I0)
    if N == 0:
        return I
    else:
        for k in range(N):
            Ik = []
            for t in range(len(x)):
                Ik.append(np.trapz(I[-1][:t+1],x[:t+1]))
            I.append(Ik)
        return np.array(I)

def Kaufmann2003Evaluate(x,params):
    a,b = params
    assert len(a) == len(b)+1
    return a[0] + np.dot(a[1:],np.exp(-np.outer(b,x)))

def Kaufmann2003DEvaluate(x,params):
    a,b = params
    assert len(a) == len(b)+1
    return np.dot(-a[1:]*b,np.exp(-np.outer(b,x)))

def Kaufmann2003Solve(N,x,y):
    """Fit N-exponential decay to a dataseries (x, y) using algorithm by
    Kaufmann, cited below.

    Parameters
    ----------
    N : float
        number of summed exponentials to fit

    x : array
        x values

    y : array
        y values

        returns a, b
        len(a) = 1 + N
        len(b) = N

        y(x) = a_0 + \sum_{i=1}^N a_i \exp ( - b_i x )

        Algorithm based on
            Kaufmann, Bernhard. “Fitting a Sum of Exponentials to Numerical Data.”
            arXiv:physics/0305019, May 6, 2003.
            http://arxiv.org/abs/physics/0305019.

    """
    x = np.array(x)
    y = np.array(y)

    # compute Ik
    Ik = ComputeIkFromY(N,x,y)

    # solve for c_i, d_i
    c_0 = np.arange(N)
    d_0 = np.arange(N+1)
    u_0 = np.array(c_0.tolist() + d_0.tolist())
    def R(u,x,N,y):
        c = u[:N]
        d = u[N:]
        return np.sum(c*Ik[1:].T,axis=1) + np.sum(d[1:]*np.array([x]*N).T**(np.arange(N)+1),axis=1) + d[0] - y

    u,ier = opt.leastsq(R,u_0,args=(x,N,y))
    c = u[:N]
    d = u[N:]

    # solve for b_i
    b = np.roots([1]+[(-1)**i*c[i] for i in range(N)])

    # solve for a_i
    a_0 = np.arange(N+1)
    def P2(a,b,d,N):
        alpha,beta = ComputeAlphaBeta(N,a,b)
        R = np.zeros(len(d))
        R[0] = d[0] - np.sum(a)
        R[1:] = d[1:]*factorial(np.arange(N)+1) - np.sum(alpha[1:]*beta[1:],axis=1)
        return R
    a = opt.fsolve(P2,a_0,(b,d,N))

    return a,b
