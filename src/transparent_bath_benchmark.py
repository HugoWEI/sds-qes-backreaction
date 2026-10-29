"""Auxiliary nongravitating smooth conformal bath appended to frozen SdS.
Not part of Easson's action; replica saddles and boundary Einstein
equations are NOT supplied. This is a free Dirac field-theory benchmark.
"""
from dataclasses import dataclass
import numpy as np
from baseline import kappa
from linearized_geometry import f0,f0_prime,f0_double_prime
from two_patch import tortoise
from free_dirac_multiinterval import _twointerval_log,_cross_ratio

@dataclass(frozen=True)
class SmoothBath:
    b:float=500.
    collar:float=50.
    mass:float=100.
    L:float=3e-6
    G:float=1.
    N:float=12.
    cutoff:float=1.

    def f(self,r):return f0(r,self.mass,self.L,self.G)
    def ygrav(self,r):
        return tortoise(r,self.b,mass=self.mass,L=self.L,G=self.G,N=self.N)
    def grav_rho(self,r):return .5*np.log(self.f(r))
    def bath_rho(self,y):
        if y<0:raise ValueError('bath y must be >=0')
        fb=self.f(self.b)
        fp=f0_prime(self.b,self.mass,self.L,self.G)
        fpp=f0_double_prime(self.b,self.mass,self.L,self.G)
        a=fp/2;d2=fb*fpp/2;z=y/self.collar
        return (.5*np.log(fb)+a*self.collar*(-np.expm1(-z))+
                (d2+a/self.collar)*self.collar**2*(1-(1+z)*np.exp(-z)))
    def collar_derivatives(self):
        fb=self.f(self.b)
        fp=f0_prime(self.b,self.mass,self.L,self.G)
        fpp=f0_double_prime(self.b,self.mass,self.L,self.G)
        return .5*np.log(fb),fp/2,fb*fpp/2
    def endpoint(self,side,t,position):
        if side=='grav':
            r=position
            if not(self.f(r)>0 and r<=self.b):
                raise ValueError('gravity endpoint exterior below b required')
            y=self.ygrav(r);rho=self.grav_rho(r)
        elif side=='bath':
            y=position;rho=self.bath_rho(y)
        else:raise ValueError('side grav/bath required')
        kb,kc=kappa(self.mass,self.L,self.G)
        u=t-y;v=t+y
        U=-np.exp(-kb*u);V=np.exp(kc*v)
        rhoUV=rho-.5*(np.log(kb*kc)-kb*u+kc*v)
        return y,U,V,rhoUV
    def _endpoints(self,points):
        z=np.array([self.endpoint(*p) for p in points],float)
        if not np.all(np.diff(z[:,0])>0):
            raise ValueError('increasing spatial y order required')
        for k in [1,2]:
            d=np.diff(z[:,k])
            if not(np.all(d>0) or np.all(d<0)):
                raise ValueError('not consistently spacelike ordered')
        return z
    def entropy_single(self,p,q):
        a,b=self._endpoints([p,q])
        return self.N/6*(np.log(abs((b[1]-a[1])*(b[2]-a[2])))+
                         a[3]+b[3]-2*np.log(self.cutoff))
    def entropy_two(self,points):
        z=self._endpoints(points)
        return self.N/6*(_twointerval_log(z[:,1])+_twointerval_log(z[:,2])+
                         z[:,3].sum()-4*np.log(self.cutoff))
    def MI(self,points):
        z=self._endpoints(points)
        xu=_cross_ratio(z[:,1]);xv=_cross_ratio(z[:,2])
        if not(0<xu<1 and 0<xv<1):raise ValueError('bad cross ratios')
        return -self.N/6*(np.log1p(-xu)+np.log1p(-xv))

if __name__=='__main__':
    m=SmoothBath()
    p=[('grav',0,320.),('grav',0,400.),('bath',0,20.),('bath',0,100.)]
    print('S(A union R)',m.entropy_two(p),'I(A:R)',m.MI(p))
    print('seam derivatives',m.collar_derivatives())
