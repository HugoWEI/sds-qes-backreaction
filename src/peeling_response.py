"""Adiabatic first-order outgoing CHIRAL interval entropy via peeling map.

Conditional state ansatz: U'(u)=exp(-integral_0^u kappa_b[M(s)]ds),
M(u)=M0-epsilon*J0*u, with FIXED u-coordinate UV cutoffs.
This omits detector redshift/counterterms, higher-order geometry and quantum
state response. It is NOT a Page curve, collector/bath, or an island saddle.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from static_patch_candidate import Model
from radiation_algebra import OutgoingAlgebra


def kappa_mass_derivative(m=Model()):
    rb,_=m.horizons()
    return -m.G*(1/rb**2+m.L)/(1-m.L*rb**2)


def alpha_slope(m=Model()):
    return -OutgoingAlgebra(m).net_flux()*kappa_mass_derivative(m)


def stationary_entropy(T,m=Model()):
    kb,_=m.kappas()
    I0=-np.expm1(-kb*T)/kb
    return m.c/6*(np.log(I0/m.uv_cutoff)+kb*T/2)


def analytic_response(T,m=Model()):
    kb,_=m.kappas()
    alpha=alpha_slope(m)
    I0=-np.expm1(-kb*T)/kb
    I2=quad(lambda x:x*x*np.exp(-kb*x),0,T,epsrel=2e-12)[0]
    return m.c/6*alpha*(T*T/4-I2/(2*I0))


def entropy_quadrature(T,epsilon=0.,m=Model()):
    kb,_=m.kappas()
    alpha=alpha_slope(m)
    K=lambda u:kb*u+epsilon*alpha*u*u/2
    I=quad(lambda u:np.exp(-K(u)),0,T,epsrel=2e-12)[0]
    return m.c/6*(np.log(I/m.uv_cutoff)+K(T)/2)


def entropy_independent_ode(T,epsilon=0.,m=Model()):
    kb,_=m.kappas()
    alpha=alpha_slope(m)
    sol=solve_ivp(lambda u,y:[-(kb+epsilon*alpha*u)*y[0],y[0]],
                  (0,T),[1.,0.],rtol=1e-12,atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    return m.c/6*(np.log(sol.y[1,-1]/m.uv_cutoff)+
                  (kb*T+epsilon*alpha*T*T/2)/2)


if __name__=="__main__":
    m=Model();T=20000.
    print("d kappa/dM, alpha",kappa_mass_derivative(m),alpha_slope(m))
    print("relative mass loss",OutgoingAlgebra(m).net_flux()*T/m.M)
    print("relative kappa change",alpha_slope(m)*T/m.kappas()[0])
    print("stationary S",stationary_entropy(T,m))
    print("analytic dS/deps",analytic_response(T,m))
    print("quadrature, independent ODE",
          entropy_quadrature(T,1,m),entropy_independent_ode(T,1,m))
