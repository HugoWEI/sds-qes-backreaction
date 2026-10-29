"""Outgoing chiral mode-interval entropy on a static observation worldtube.

Classical instantaneous SdS + conformal U=-exp(-kappa_b*u) vacuum.
This is an outgoing mode subalgebra entropy, NOT a dynamical collected
bath entropy, replica island, or physical fine-grained Page curve.
"""
from dataclasses import dataclass
import numpy as np
from static_patch_candidate import Model


@dataclass(frozen=True)
class OutgoingAlgebra:
    model: Model = Model()

    def entropy(self,delta_t):
        if delta_t<=0:
            raise ValueError("interval length must be positive")
        m=self.model
        kb,_=m.kappas()
        x=kb*delta_t/2
        logsinh=x+np.log1p(-np.exp(-2*x))-np.log(2)
        return m.c/6*(np.log(2*np.sqrt(m.f(m.b))/
                             (kb*m.uv_cutoff))+logsinh)

    def rate(self,delta_t):
        if delta_t<=0:
            raise ValueError("interval length must be positive")
        kb,_=self.model.kappas()
        return self.model.c*kb/(12*np.tanh(kb*delta_t/2))

    def asymptotic_rate(self):
        kb,_=self.model.kappas()
        return self.model.c*kb/12

    def net_flux(self):
        m=self.model;kb,kc=m.kappas()
        return m.c*(kb*kb-kc*kc)/(48*np.pi)

    def thermodynamic_proxy_rate(self):
        kb,kc=self.model.kappas()
        return self.net_flux()*2*np.pi*(1/kc-1/kb)


if __name__=="__main__":
    alg=OutgoingAlgebra()
    print("outgoing chiral late-time entropy growth",alg.asymptotic_rate())
    print("Easson thermodynamic proxy entropy growth",alg.thermodynamic_proxy_rate())
    for t in (100.,500.,1000.,5000.):
        print("duration, regulated S_out, rate",t,alg.entropy(t),alg.rate(t))
