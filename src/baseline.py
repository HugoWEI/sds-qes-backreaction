"""Easson thermodynamic Page proxy, NOT a microscopic QES computation.
Units: c=hbar=1. Stay away from Nariai and the quantum endpoint.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


def horizons(m, L=3e-6, G=1.0):
    assert 0 < 9*(G*m)**2*L < 1
    f=lambda r: 1-2*G*m/r-L*r*r/3
    turn=1/np.sqrt(L)
    return brentq(f,G*m,turn), brentq(f,turn,np.sqrt(3/L))


def kappa(m, L=3e-6, G=1.0):
    rb,rc=horizons(m,L,G)
    return (1-L*rb*rb)/(2*rb), (L*rc*rc-1)/(2*rc)


def flux(m,L=3e-6,G=1.0,N=12):
    kb,kc=kappa(m,L,G)
    return N/(48*np.pi)*(kb*kb-kc*kc)


def entropies(m,L=3e-6,G=1.0):
    rb,rc=horizons(m,L,G)
    return np.pi*rb*rb/G, np.pi*rc*rc/G


def page_mass_proxy(m0,L=3e-6,G=1.0):
    rb,rc=horizons(m0,L,G)
    rcp=np.hypot(rb,rc)
    return rcp*(1-L*rcp**2/3)/(2*G)


def page_time_proxy(m0,L=3e-6,G=1.0,N=12,eps=0.):
    mp=page_mass_proxy(m0,L,G)
    return quad(lambda m: 1/(flux(m,L,G,N)*(1+eps*m/m0)),
                mp,m0,epsabs=1e-9,epsrel=1e-11)[0]


def page_time_derivative(m0,L=3e-6,G=1.0,N=12):
    mp=page_mass_proxy(m0,L,G)
    return -quad(lambda m: (m/m0)/flux(m,L,G,N),
                 mp,m0,epsabs=1e-9,epsrel=1e-11)[0]


def check_ode(m0=100.,L=3e-6,G=1.0,N=12,eps=0.01):
    mp=page_mass_proxy(m0,L,G)
    target=page_time_proxy(m0,L,G,N,eps)
    def event(t,y):
        return y[0]-mp
    event.terminal=True
    event.direction=-1
    sol=solve_ivp(lambda t,y: [-flux(y[0],L,G,N)*(1+eps*y[0]/m0)],
                  (0,2*target),[m0],events=event,rtol=1e-11,atol=1e-12)
    return target,sol.t_events[0][0]


if __name__=="__main__":
    M0=100.
    rb,rc=horizons(M0); mp=page_mass_proxy(M0)
    s0=sum(entropies(M0)); sb,sc=entropies(mp)
    integral=quad(lambda m: 2*np.pi*(1/kappa(m)[1]-1/kappa(m)[0]),
                  mp,M0)[0]
    h=1e-4
    fd=(page_time_proxy(M0,eps=h)-page_time_proxy(M0,eps=-h))/(2*h)
    tq,to=check_ode()
    print("initial horizons",rb,rc)
    print("proxy Page mass/time",mp,page_time_proxy(M0))
    print("entropy integral/DeltaStotal/S_b",integral,sb+sc-s0,sb)
    print("dt/depsilon analytic/finite difference",page_time_derivative(M0),fd)
    print("perturbed time quadrature/independent ODE",tq,to)
