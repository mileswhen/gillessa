import numpy as np
from numpy.random import exponential as rexp


class SSA:
    def __init__(self, S, props, x0):
        '''
        STOCHASTIC SIMULATION ALGORITHM
        doi: 10.1021/j100540a008
        ---
        S : nd_array
            transpose of stoich matrix of reaction species coefficients
        props : dict 
            reaction propensity functions
        x0 : dict
            init number of species 
        '''
        assert(S.shape[0] == len(props))
        
        self.S = S           
        self.props = props
        self.X = x0
        
    def __w(self, m):
        '''convenience method for getting propensity w_k'''
        wk = self.props[m](self.X)

        return wk
    
    def __getk(self, w0):
        u = np.random.uniform()
        p = 0
        for m in range(1, len(self.props)+1):
            pm = self.__w(m) / w0 
            if p < u <= p+pm:
                k = m
            p += pm
        
        return k
        
    def __update_state(self, k):
        '''update state'''
        for species, coeff in zip(self.X, self.S[k-1]):
            self.X[species] += coeff
            
    def sim(self, T_end):
        states = [list(self.X.values())]
        times = [0]
        t = 0
        
        while t <= T_end:
            # compute w0
            w0 = 0
            for k in range(1, len(self.props)+1):
                w0 += self.__w(k)
                
            # sample T0
            T0 = rexp(1/w0)
            if t + T0 > T_end:
                break
            else:
                t += T0
            
            # choose reaction and update state
            k = self.__getk(w0)
            self.__update_state(k)
            
            times.append(t)
            states.append(list(self.X.values()))
            
        return times, states


class FRM:
    def __init__(self, S, props, x0):
        '''
        FIRST REACTION METHOD
        doi: 10.1016/0021-9991(76)90041-3
        ---
        S : nd_array
            transpose of stoich matrix of reaction species coefficients
        props : dict 
            reaction propensity functions
        x0 : dict
            init number of species 
        '''
        assert(S.shape[0] == len(props))
        
        self.S = S           # (row=reactions, col=species)
        self.props = props
        self.X = x0
        
    def __w(self, m):
        '''convenience method for reaction propensity w_k'''
        wk = self.props[m](self.X)

        return wk
    
    def __sampleExp(self, k):
        '''samples time_k ~ exponential'''
        wk = self.__w(k)
        if wk != 0:
            tk = rexp(1/wk)
        else:
            tk = np.inf

        return tk
    
    def __get_tk(self):
        '''choose reaction'''
        t_ks = {k: self.__sampleExp(k) for k in range(1, len(self.props)+1)}
        k = min(t_ks, key=t_ks.get)
        
        return t_ks[k], k
    
    def __update_state(self, k):
        '''update state'''
        for species, coeff in zip(self.X, self.S[k-1]):
            self.X[species] += coeff
    
    def sim(self, T_end):
        states = [list(self.X.values())]
        times = [0]
        t = 0
        
        while t <= T_end:
            tk, k = self.__get_tk()
            if t + tk > T_end:
                break
                
            t += tk
            self.__update_state(k)
            
            times.append(t)
            states.append(list(self.X.values()))

        return times, states 


class MNRM:
    def __init__(self, S, props, x0):
        '''
        MODIFIED NEXT REACTION METHOD
        doi: 10.1063/1.2799998
        ---
        S : nd_array
            transpose of stoich matrix of reaction species coefficients
        props : dict 
            reaction propensity functions
        x0 : dict
            init number of species 
        '''
        assert(S.shape[0] == len(props))
        
        self.S = S           # (row=reactions, col=species)
        self.props = props
        self.X = x0
        
    def __w(self, m):
        '''convenience method for getting propensity w_k'''
        wk = self.props[m](self.X)

        return wk
    
    def __get_delta(self, P, T):
        '''convenience method for getting delta_k'''
        delta_k = {}
        for k in range(1, len(P)+1):
            wk = self.__w(k)
            if wk != 0:
                delta_k[k] = (P[k] - T[k]) / wk
            else:
                delta_k[k] = np.inf
        
        k = min(delta_k, key=delta_k.get)
        
        return delta_k[k], k
    
    def __update_state(self, k):
        '''update state'''
        for species, coeff in zip(self.X, self.S[k-1]):
            self.X[species] += coeff
            
    def sim(self, T_end):
        states = [list(self.X.values())]
        times = [0]
        P = {k: rexp(1) for k in range(1, len(self.props)+1)}
        T = {k: 0 for k in range(1, len(self.props)+1)}
        t = 0

        while t < T_end:
            # delta = min_k{delta_k}
            delta, k = self.__get_delta(P, T)
            if t + delta > T_end:
                break

            # update T
            for m in range(1, len(self.props)+1):
                T[m] += self.__w(m) * delta

            # update Pk
            P[k] += rexp(1)
            
            # update state and time
            self.__update_state(k)
            t += delta
            times.append(t)
            states.append(list(self.X.values()))
            
        return times, states


class TAULP:
    def __init__(self, S, props, x0, tau):
        '''
        TAUL LEAP METHOD
        doi: 10.1063/1.1378322
        ---
        S : nd_array
            transpose of stoich matrix of reaction species coefficients
        props : dict 
            reaction propensity functions
        x0 : dict
            init number of species 
        '''
        assert(S.shape[0] == len(props))
        
        self.S = S          
        self.props = props
        self.X = x0
        self.tau = tau
        
    def __w(self, k):
        '''convenience method for getting propensity w_k'''
        wk = self.props[k](self.X)

        return wk
    
    def __get_Nk(self, k):
        lamb = self.__w(k) * self.tau
        if lamb >= 0:
            Nk = np.random.poisson(lamb)
        else:
            Nk = 0
        return Nk
    
    def __update_state(self, Nk, k):
        '''update state'''
        if Nk != 0:
            for species, coeff in zip(self.X, self.S[k-1]):
                self.X[species] += coeff * Nk
        
    def sim(self, T_end):
        states = [list(self.X.values())]
        times = [0]
        t = 0
        
        while t <= T_end:
            if t + self.tau > T_end:
                break
            else:
                t += self.tau
            
            # update state for all reactions
            Nks = {k: self.__get_Nk(k) for k in range(1, len(self.props)+1)}
            for k in range(1, len(self.props)+1):
                self.__update_state(Nks[k], k)
                
            # remove negatives to prevent numerical issues
            for species, n in self.X.items():
                if n < 0:
                    self.X[species] = 0

            times.append(t)
            states.append(list(self.X.values()))            
            
        return times, states
