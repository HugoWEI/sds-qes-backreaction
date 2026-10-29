"""Conditional 2D CFT interval entropy and linearized response in the overlap.

Assumptions: Gaussian conformal vacuum in U=-exp(-kb*u), V=-exp(-kc*v),
physical UV cutoffs fixed, two endpoints in static region, separated
spacelike. Does NOT implement a replica island or Easson Page entropy proxy.
"""
import numpy as np
from baseline import kappa, flux
from linearized_geometry import f0, radial_solution
from two_patch import tortoise, outgoing_radial_solution, coordinate_correction


def background_uv(v,r,r_reference=400.,**kw):
    return v-2*tortoise(r,r_reference,**kw),v


def linearized_rho(v,r, r_reference=400., **kw):
    """rho = (1/2)log(f0) + eps*rho1 in coordinate (u,v)."""
    sigp,Fr=radial_solution(r,r_reference,**kw)
    sigm,_=outgoing_radial_solution(r,r_reference,**kw)
    G=kw.get('G',1.)
    J=flux(kw.get('mass',100.),kw.get('L',3e-6),G,kw.get('N',12))
    F1=Fr+2*G*J*v/r
    f=f0(r,kw.get('mass',100.),kw.get('L',3e-6),G)
    return 0.5*np.log(f),0.5*(F1/f+sigp+sigm)


def interval_entropy(v1,r1,v2,r2,cutoff1=1.,cutoff2=1.,
                     epsilon=0.,r_reference=400.,**kw):
    """Renormalized interval entropy, with chosen cutoff lengths; not Page."""
    N=kw.get('N',12)
    kb,kc=kappa(kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
    u1,w1=background_uv(v1,r1,r_reference,**kw)
    u2,w2=background_uv(v2,r2,r_reference,**kw)
    u1+=epsilon*coordinate_correction(v1,r1,r_reference,**kw)
    u2+=epsilon*coordinate_correction(v2,r2,r_reference,**kw)
    du=u1-u2;dv=w1-w2
    if du*dv>=0:
        raise ValueError('Endpoints must be spacelike in (u,v)')
    rho01,rho11=linearized_rho(v1,r1,r_reference,**kw)
    rho02,rho12=linearized_rho(v2,r2,r_reference,**kw)
    uv_product=np.abs((2/kb)*np.sinh(kb*du/2)*(2/kc)*np.sinh(kc*dv/2))
    return N/6*(np.log(uv_product/(cutoff1*cutoff2))+
                rho01+rho02+epsilon*(rho11+rho12))


def interval_entropy_response(v1,r1,v2,r2,r_reference=400.,**kw):
    """Analytic epsilon derivative of the stipulated two-point entropy."""
    N=kw.get('N',12);kb,_=kappa(kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
    u1,_=background_uv(v1,r1,r_reference,**kw)
    u2,_=background_uv(v2,r2,r_reference,**kw)
    chi1=coordinate_correction(v1,r1,r_reference,**kw)
    chi2=coordinate_correction(v2,r2,r_reference,**kw)
    _,rho11=linearized_rho(v1,r1,r_reference,**kw)
    _,rho12=linearized_rho(v2,r2,r_reference,**kw)
    return N/6*(kb/2 / np.tanh(kb*(u1-u2)/2)*(chi1-chi2)+rho11+rho12)


if __name__=='__main__':
    from two_patch import tortoise
    dv=tortoise(500.)
    print('stationary entropy for equal-static-times:')
    for t in [0.,100.,1000.]:
        print(t,interval_entropy(t,400.,t+dv,500.))
    s1=interval_entropy_response(0.,400.,dv,500.)
    print('linearized interval entropy response:',s1)
    for h in [0.5,0.1,0.01]:
        fd=(interval_entropy(0.,400.,dv,500.,epsilon=h)-interval_entropy(0.,400.,dv,500.,epsilon=-h))/(2*h)
        print('h',h,'central-difference',fd,'residual',fd-s1)
