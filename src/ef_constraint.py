"""Local O(NG) radial constraint for generalized EF metric.

This is a frozen-stress diagnostic, not a global UdS QES or complete metric.
Two-dimensional Polyakov stress and G conventions are explicit.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from baseline import horizons,kappa


def blackening(r,m=100.,L=3e-6,G=1.):
    return 1-2*G*m/r-L*r*r/3


def static_chiral_tuu(r,m=100.,L=3e-6,G=1.,N=12):
    f=blackening(r,m,L,G)
    fp=2*G*m/r**2-2*L*r/3
    fpp=-4*G*m/r**3-2*L/3
    kb,_=kappa(m,L,G)
    tu=kb*kb/4
    return -N/(12*np.pi)*(fp*fp/16-f*fpp/8-tu)


def ef_trr(r,m=100.,L=3e-6,G=1.,N=12):
    f=blackening(r,m,L,G)
    return 4*static_chiral_tuu(r,m,L,G,N)/f**2


def sigma_prime(r,m=100.,L=3e-6,G=1.,N=12):
    # E_rr = 2r sigma' = 2G T_rr
    return G*ef_trr(r,m,L,G,N)/r


def sigma_difference(r0,r1,m=100.,L=3e-6,G=1.,N=12):
    return quad(lambda r:sigma_prime(r,m,L,G,N),r0,r1,
                epsabs=1e-12,epsrel=1e-11)[0]


def sigma_difference_ode(r0,r1,m=100.,L=3e-6,G=1.,N=12):
    sol=solve_ivp(lambda r,y:[sigma_prime(r,m,L,G,N)],
                  (r0,r1),[0.],rtol=1e-11,atol=1e-13)
    assert sol.success
    return sol.y[0,-1]


if __name__=="__main__":
    rb,rc=horizons(100.)
    r0=(rb+rc)/2; r1=500.
    print("initial horizons",rb,rc)
    print("Tuu, Trr at r=500",static_chiral_tuu(r1),ef_trr(r1))
    print("delta sigma via quadrature / ODE",
          sigma_difference(r0,r1),sigma_difference_ode(r0,r1))
