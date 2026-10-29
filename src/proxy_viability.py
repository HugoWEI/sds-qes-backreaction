"""Necessary instantaneous semiclassical diagnostics along Easson's Page PROXY.

This is NOT an island, a fine-grained Page curve, or a bath validity test.
Use only sub-Nariai neutral SdS, with dimensionless x=sqrt(Lambda)*r_b,
q=N*G*Lambda. Flux normalization is Easson's declared model convention.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq


def cosmological_y(x):
    if not 0 < x < 1:
        raise ValueError('0 < sqrt(Lambda) r_b < 1 required')
    return (np.sqrt(12-3*x*x)-x)/2


def mass_mu(x):
    """Dimensionless mass mu=G M sqrt(Lambda)."""
    return x*(1-x*x/3)/2


def kappas_scaled(x):
    y=cosmological_y(x)
    return (1-x*x)/(2*x),(y*y-1)/(2*y)


def proxy_page(x):
    """Exact Easson proxy crossing, NOT a QES Page time.

    Entropy identity: y_page^2=x_initial^2 + y_initial^2.
    With x^2+x*y+y^2=3, p=x*y, hence mu_page=p*sqrt(3-p)/6.
    """
    y=cosmological_y(x)
    p=x*y
    mu_page=p*np.sqrt(3-p)/6
    xp=brentq(lambda v:mass_mu(v)-mu_page,1e-15,x,xtol=1e-14)
    return dict(x_page=xp,mu_page=mu_page,initial_mu=mass_mu(x),
                mass_ratio=mu_page/mass_mu(x),
                initial_p=p, initial_y=y)


def adiabatic_coefficient(x):
    """A(x): |dot(kappa_b)|/kappa_b^2 = q*A(x)/(48*pi).

    Assumes dot M=-N/(48*pi)*(kappa_b^2-kappa_c^2),
    d kappa_b/dM= -G*(1/r_b^2+Lambda)/(1-Lambda*r_b^2).
    """
    if x>1-1e-5:
        # Analytic regular Nariai expansion; avoids 0/0 cancellation.
        return 8/3+56*(x-1)**2/9
    y=cosmological_y(x)
    kb,kc=kappas_scaled(x)
    return (1+x*x)/(x*x*(1-x*x))*(1-(kc/kb)**2)


def cosmological_adiabatic_coefficient(x):
    """B(x): |dot(kappa_c)|/kappa_c^2 = q*B(x)/(48*pi).

    This is the second-horizon analogue of adiabatic_coefficient.
    """
    if x>1-1e-5:
        return 8/3  # Nariai limit; corrections vanish as x->1.
    y=cosmological_y(x)
    kb,kc=kappas_scaled(x)
    return (1+y*y)/(y*y*(y*y-1))*((kb/kc)**2-1)


@dataclass(frozen=True)
class ProxyDiagnostic:
    initial_x: float
    q: float = 1e-6  # q=N G Lambda, not independently a 4D bath coupling.

    def evaluate(self,grid=501):
        if self.q<=0:raise ValueError('q>0 required')
        page=proxy_page(self.initial_x)
        xp=page['x_page']
        xs=np.linspace(xp,self.initial_x,grid)
        a=np.array([adiabatic_coefficient(x) for x in xs])
        b=np.array([cosmological_adiabatic_coefficient(x) for x in xs])
        coupling=self.q/(xs*xs) # cG/r_b^2, ignoring O(1) constants
        return dict(**page,initial_x=self.initial_x,q=self.q,
                    max_adiabatic=float(self.q*np.max(np.maximum(a,b))/(48*np.pi)),
                    max_bh_adiabatic=float(self.q*np.max(a)/(48*np.pi)),
                    max_cosmological_adiabatic=float(self.q*np.max(b)/(48*np.pi)),
                    max_coupling=float(np.max(coupling)),
                    eps_at_page=float(self.q*adiabatic_coefficient(xp)/(48*np.pi)),
                    critical_q_for_one_percent=float(.01*48*np.pi/np.max(np.maximum(a,b))))


if __name__=='__main__':
    print('Dimensionless necessary validity diagnostics along THERMODYNAMIC proxy')
    print('Nariai A limit = 8/3; small-x x^2*A -> 1; neither is a QES theorem')
    for x in [.05,.2,.5,.8,.99]:
        d=ProxyDiagnostic(x).evaluate()
        print('x0=',x,'Mp/M0=',format(d['mass_ratio'],'.8f'),
              'xPage=',format(d['x_page'],'.8f'),
              'max adiab=',format(d['max_adiabatic'],'.8g'),
              'max cG/rb^2=',format(d['max_coupling'],'.8g'),
              'qcrit(eps<.01)=',format(d['critical_q_for_one_percent'],'.8g'))
