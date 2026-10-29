"""Two-null stationary point of a UV-invalid exterior TRIAL functional.
This does NOT produce a genuine island/replica saddle or a Page time.
"""
import numpy as np
from scipy.optimize import root
from scipy.integrate import quad
from transparent_bath_benchmark import SmoothBath
from baseline import kappa
from linearized_geometry import f0_prime

def sample_points(a,t,m=SmoothBath(),fixed=400.,R1=20.,R2=100.):
    return [('grav',t,a),('grav',0.,fixed),('bath',0.,R1),('bath',0.,R2)]

def trial_entropy(a,t,m=SmoothBath(),**kw):
    return np.pi*a*a/m.G+m.entropy_two(sample_points(a,t,m,**kw))

def analytic_gradient(a,t,m=SmoothBath(),**kw):
    z=m._endpoints(sample_points(a,t,m,**kw))
    U=z[:,1];V=z[:,2]
    A=lambda v:1/(v[0]-v[1])+1/(v[0]-v[3])-1/(v[0]-v[2])
    kb,kc=kappa(m.mass,m.L,m.G)
    p=kb*U[0]*A(U);q=kc*V[0]*A(V)
    f=m.f(a);fp=f0_prime(a,m.mass,m.L,m.G)
    return np.array([2*np.pi*a/m.G+m.N/(6*f)*(p+q+(fp-kb-kc)/2),
                     m.N/6*(-p+q+(kb-kc)/2)])

def find_root(m=SmoothBath(),fixed=400.,**kw):
    gap0=m.N*m.G/(6*np.pi*fixed)
    def fun(z):
        try:return analytic_gradient(z[0],z[1],m,fixed=fixed,**kw)
        except (ValueError,OverflowError,FloatingPointError):
            return np.array([1e50,1e50])
    sol=root(fun,np.array([fixed-gap0,0.]),method='hybr',tol=1e-10)
    if not sol.success and np.linalg.norm(fun(sol.x))>1e-4:
        raise RuntimeError(sol.message)
    return sol.x[0],sol.x[1],fun(sol.x)

def proper_gap(a,fixed=400.,m=SmoothBath()):
    return quad(lambda r:1/np.sqrt(m.f(r)),a,fixed)[0]

if __name__=='__main__':
    m=SmoothBath();a,t,g=find_root(m)
    print('formal a,t',a,t,'residual',g)
    print('radial gap',400-a,'proper gap',proper_gap(a,m=m),
          'cutoff',m.cutoff,'near-UV estimate',m.N*m.G/(6*np.pi*400))
