"""Two-null extremization stress-test on frozen SdS, NOT a physical island.

Assumptions: S=pi*a^2/G+S_CFT([island_endpoint, fixed_endpoint_b])
for one exterior spacelike interval, U=-exp(-kb*u), V=-exp(-kc*v).
A physical bath/collector, purification and replica channel ARE NOT provided.
The only stationary point in the reference parameter scan is UV-adjacent
to the fixed endpoint, not near the black hole horizon.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Model:
    G: float = 1.
    L: float = 3e-6
    M: float = 100.
    c: float = 12.
    b: float = 500.
    uv_cutoff: float = 1.

    def f(self, r):
        return 1-2*self.G*self.M/r-self.L*r*r/3

    def fp(self, r):
        return 2*self.G*self.M/r**2-2*self.L*r/3

    def horizons(self):
        turn=1/np.sqrt(self.L)
        return (brentq(self.f,self.G*self.M,turn),
                brentq(self.f,turn,np.sqrt(3/self.L)))

    def kappas(self):
        rb,rc=self.horizons()
        return self.fp(rb)/2,-self.fp(rc)/2

    def rstar(self,r):
        rb,rc=self.horizons()
        ro=-rb-rc
        return sum(np.log(abs((r-h)/(self.b-h)))/self.fp(h)
                   for h in (rb,rc,ro))

    def effective_D(self,a):
        assert self.horizons()[0] < a < self.b < self.horizons()[1]
        return -self.rstar(a)

    def time_gradient(self,a,dt):
        kb,kc=self.kappas();D=self.effective_D(a)
        x=dt+D; y=D-dt
        if min(x,y)<=0:
            return np.nan
        return self.c/12*(kb/np.tanh(kb*x/2)-kc/np.tanh(kc*y/2))

    def time_stationary(self,a):
        D=self.effective_D(a)
        tol=max(1e-14,1e-12*D)
        return brentq(lambda dt:self.time_gradient(a,dt),
                      -D+tol,D-tol,xtol=1e-12)

    def radial_gradient(self,a,dt):
        kb,kc=self.kappas();D=self.effective_D(a)
        x=dt+D; y=D-dt
        if min(x,y)<=0:
            return np.nan
        c1=kb/np.tanh(kb*x/2)
        c2=kc/np.tanh(kc*y/2)
        return 2*np.pi*a/self.G+self.c/(12*self.f(a))*(self.fp(a)-c1-c2)

    def generalized_entropy(self,a,dt):
        kb,kc=self.kappas();D=self.effective_D(a)
        x=dt+D;y=D-dt
        if min(x,y)<=0:
            return np.inf
        value=abs(4*np.sinh(kb*x/2)*np.sinh(kc*y/2)/(kb*kc))
        return (np.pi*a*a/self.G+self.c/6*np.log(
                value*np.sqrt(self.f(a)*self.f(self.b))/self.uv_cutoff**2))

    def reduced_gradient(self,a):
        return self.radial_gradient(a,self.time_stationary(a))

    def find_near_endpoint_root(self,search_window=1.):
        gap=brentq(lambda d:self.reduced_gradient(self.b-d),
                   1e-6,search_window,xtol=1e-13)
        a=self.b-gap
        return a,self.time_stationary(a)

    def uv_scale(self,a):
        return (self.b-a)/np.sqrt(self.f(self.b))

    def hessian(self,a,dt):
        ha=1e-4*(self.b-a)
        ht=ha/self.f(self.b)
        S=self.generalized_entropy
        Haa=(S(a+ha,dt)-2*S(a,dt)+S(a-ha,dt))/ha**2
        Htt=(S(a,dt+ht)-2*S(a,dt)+S(a,dt-ht))/ht**2
        Hat=(S(a+ha,dt+ht)-S(a+ha,dt-ht)
             -S(a-ha,dt+ht)+S(a-ha,dt-ht))/(4*ha*ht)
        return np.array([[Haa,Hat],[Hat,Htt]])


if __name__=="__main__":
    m=Model()
    a,dt=m.find_near_endpoint_root()
    print("formal two-null endpoint:",a,dt)
    print("gap / leading UV estimate",m.b-a,m.c*m.G/(6*np.pi*m.b))
    print("proper separation / UV cutoff",m.uv_scale(a)/m.uv_cutoff)
    print("two gradients",m.time_gradient(a,dt),m.radial_gradient(a,dt))
    print("Hessian eigenvalues",np.linalg.eigvalsh(m.hessian(a,dt)))
