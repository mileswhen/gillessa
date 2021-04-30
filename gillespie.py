import numpy as np

# not an object-oriented approach, but who has time for that anyways...

def tauTime(avec):
    '''returns tau given r1 in (0,1)'''
    r1 = np.random.uniform()
    a0 = sum(avec)
    # tau = (1/a0)*np.log(1/r1)
    tau = np.random.exponential(1/a0)

    return tau


def get_alphas(pvec, kvec):
    '''Returns probability of reaction events: alpha_i
    pvec : [a_t, f_t, fa_t, faf_t]
    kvec : [k1, k2, k3, k4]
    '''
    a_t, f_t, fa_t, faf_t = pvec
    k1, k2, k3, k4 = kvec

    # propensity functions
    a1 = f_t * a_t * k1
    a2 = fa_t * k2
    a3 = fa_t * f_t * k3
    a4 = faf_t * k4

    return [a1, a2, a3, a4]


def chooseR(pvec, avec)
    '''Chooses reaction'''
    r2 = np.random.uniform()
    a_t, f_t, fa_t, faf_t = pvec
    a1, a2, a3, a4 = avec
    a0 = sum(avec)

    # a1 conditions
    if r2 < a1/a0:
        a_t -= 1
        f_t -= 1
        fa_t += 1
    
    # a2 conditions
    if a1/a0 <= r2 < (a1+a2)/a0:
        a_t += 1
        f_t += 1
        fa_t -= 1
    
    # a3 conditions:
    if (a1+a2)/a0 <= r2 < (a1+a2+a3)/a0:
        f_t -= 1
        fa_t -= 1
        faf_t += 1
    
    # a4 conditions:
    if (a1+a2+a3)/a0 <= r2 < (a1+a2+a3+a4)/a0:
        f_t += 1
        fa_t += 1
        faf_t -= 1
    
    return [a_t, f_t, fa_t, faf_t]


def sim(pvec_initial, kvec, maxstep):
    '''Returns list of lists for particle number evolution: [At, Ft, FAt, FAFt]'''
    plist = [pvec_initial]
    tlist = [0]

    for step in range(1, maxstep):
        '''simulate'''
        # vector of particle numbers at step: [a_t, f_t, fa_t, faf_t]
        pvec = plist[step-1]

        # get vector of reaction probabilities avec
        avec = get_alphas(pvec=pvec, kvec=kvec)

        # add tau time to list
        tau = tauTime(avec)
        t = tlist[step-1]
        tlist.append(t+tau)

        # choose reaction and update pvec
        pvec_update = chooseR(pvec, avec)
        plist.append(pvec_update)
        
    return plist, tlist

