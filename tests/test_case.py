import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest
from models import SSA, FRM, MNRM, TAULP


def parameterize():
    # parametrize simulation:
    S = np.array([[-1, 1, 0, 0],
                  [1, -1, 0, 0],
                  [-2, 0, 1, 0],
                  [2, 0, -1, 0],
                  [1, 0, 0, -1], 
                  [-1, 0, 0, 1]])
    props = {
        1: lambda X: X['X'],
        2: lambda X: X['Y'],
        3: lambda X: X['X']*(X['X']-1),
        4: lambda X: 2*X['Z'],
        5: lambda X: 0.1*X['X']*X['W'],
        6: lambda X: 0.05*X['X']*(X['X']-1)
    }
    x0 = {'X': 20, 'Y': 0, 'Z': 0, 'W': 10}
    return S, props, x0


def get_moments(states):
    mu = np.mean(states, axis=0)
    var = np.var(states, axis=0)  
    CV = np.sqrt(var) / mu
    return mu, var, CV


# run tests
def test_ssa():
    # simulate
    states = []
    S, props, x0 = parameterize()
    for _ in range(10000):
        ssa = SSA(S, props, 
                  {'X': 20, 'Y': 0, 'Z': 0, 'W': 10})
        times, x = ssa.sim(1)
        states.append(x[-1])

    # mu, std, CV
    mu, var, CV = get_moments(states)
    assert(4 < mu[0] < 4.15)
    assert(2.6 < mu[1] < 2.8)
    assert(7.9 < mu[2] < 8.1)
    assert(7.15 < mu[3] < 7.35)


def test_frm():
    # simulate
    states = []
    for _ in range(10000):
        S, props, x0 = parameterize()
        frm = FRM(S, props, x0)
        times, x = frm.sim(1)
        states.append(x[-1])

    # mu, std, CV
    mu, var, CV = get_moments(states)
    assert(4 < mu[0] < 4.15)
    assert(2.6 < mu[1] < 2.8)
    assert(7.9 < mu[2] < 8.1)
    assert(7.15 < mu[3] < 7.35)


def test_mnrm():
    # simulate
    states = []
    for _ in range(10000):
        S, props, x0 = parameterize()
        mnrm = MNRM(S, props, x0)
        times, x = mnrm.sim(1)
        states.append(x[-1])

    # mu, std, CV
    mu, var, CV = get_moments(states)
    assert(4 < mu[0] < 4.15)
    assert(2.6 < mu[1] < 2.8)
    assert(7.9 < mu[2] < 8.1)
    assert(7.15 < mu[3] < 7.35)


def test_taulp():
    # simulate
    states = []
    for _ in range(10000):
        S, props, x0 = parameterize()
        taulp = TAULP(S, props, x0, tau=1e-3)
        times, x = taulp.sim(1)
        states.append(x[-1])

    # mu, std, CV
    mu, var, CV = get_moments(states)
    assert(4 < mu[0] < 4.15)
    assert(2.6 < mu[1] < 2.8)
    assert(7.9 < mu[2] < 8.1)
    assert(7.15 < mu[3] < 7.35)
