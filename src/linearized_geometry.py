"""Local first-order Einstein-dilaton-Polyakov metric, NOT a global QES solution.
Only use between distinct SdS horizons in an ingoing EF chart.
The state uses chiral t_u=kappa_b^2/4, t_v=kappa_c^2/4.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from baseline import kappa, flux


def f0(r,mass=100.,L=3e-6,G=1.):
    return 1-2*G*mass/r-L*r*r/3


def f0_prime(r,mass=100.,L=3e-6,G=1.):
    return 2*G*mass/r**2-2*L*r/3


def f0_double_prime(r,mass=100.,L=3e-6,G=1.):
    return -4*G*mass/r**3-2*L/3


def stress(r,mass=100.,L=3e-6,G=1.,N=12):
    """2D conformal and ingoing-EF tensor components, on classical f0."""
    f=f0(r,mass,L,G)
    if f<=0:
        raise ValueError("This patch requires r_b<r<r_c")
    fp=f0_prime(r,mass,L,G)
    fpp=f0_double_prime(r,mass,L,G)
    kb,kc=kappa(mass,L,G)
    Tuu=-N/(192*np.pi)*(fp**2-2*f*fpp-4*kb**2)
    Tvv=-N/(192*np.pi)*(fp**2-2*f*fpp-4*kc**2)
    Tuv=N*f*fpp/(96*np.pi)
    return dict(Tuu=Tuu,Tvv=Tvv,Tuv=Tuv,
                Trr=4*Tuu/f**2,
                Tvr=-2*(Tuu+Tuv)/f,
                Tvv_ef=Tuu+2*Tuv+Tvv)


def sigma_source(r,**kw):
    return kw.get("G",1.)*stress(r,**kw)["Trr"]/r


def radial_F1_source(r,**kw):
    return 2*kw.get("G",1.)*stress(r,**kw)["Tvr"]


def radial_solution(r,r_reference=400.,**kw):
    """Quadratures, with sigma1=F1=0 at r_reference and reference time."""
    sigma1=quad(lambda x:sigma_source(x,**kw),r_reference,r,
                epsabs=1e-13,epsrel=1e-11)[0]
    F1=quad(lambda x:radial_F1_source(x,**kw),r_reference,r,
            epsabs=1e-13,epsrel=1e-11)[0]/r
    return sigma1,F1


def radial_solution_ode(r,r_reference=400.,**kw):
    """Independent numerical integration of the two first-order ODEs."""
    sol=solve_ivp(lambda x,y:[sigma_source(x,**kw),
                 (radial_F1_source(x,**kw)-y[1])/x],
                 (r_reference,r),[0.,0.],rtol=1e-12,atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    return tuple(sol.y[:,-1])


def time_dependent_F1(v,r,v_reference=0.,r_reference=400.,**kw):
    """Linear secular perturbation; |delta M|<<M is required."""
    J=flux(kw.get("mass",100.),kw.get("L",3e-6),
           kw.get("G",1.),kw.get("N",12))
    return radial_solution(r,r_reference,**kw)[1]+2*kw.get("G",1.)*J*(v-v_reference)/r


def local_metric_check(r=500.,r_reference=400.,**kw):
    """Exact O(NG) EF metric-equation residual check from source stresses."""
    G=kw.get("G",1.); m=kw.get("mass",100.)
    L=kw.get("L",3e-6);N=kw.get("N",12)
    f=f0(r,m,L,G);t=stress(r,**kw)
    sigma1,F1=radial_solution(r,r_reference,**kw)
    Err=2*r*sigma_source(r,**kw)
    Evr=radial_F1_source(r,**kw)
    Evv=-f*Evr-2*G*flux(m,L,G,N)
    return dict(E_rr=Err,T_rr=2*G*t["Trr"],
                E_vr=Evr,T_vr=2*G*t["Tvr"],
                E_vv=Evv,T_vv=2*G*t["Tvv_ef"],
                sigma1=sigma1,F1_radial=F1)


if __name__=="__main__":
    for r in [400.,500.,600.]:
        print(r,local_metric_check(r))
        print("quadrature/ODE",radial_solution(r),radial_solution_ode(r))
