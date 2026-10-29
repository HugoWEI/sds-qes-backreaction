"""Dual horizon-adapted EF charts on an instantaneous SdS background.

The perturbation is evaluated on the *same* fixed zeroth-order Polyakov state.
The coordinate matching is valid only to first order in the anomaly stress,
away from horizon boundary layers. This is not a QES computation.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from baseline import flux, horizons
from linearized_geometry import f0, stress, radial_solution


def outgoing_stress(r, **kw):
    """Tensor transformation with v=u+2*rstar(r), u fixed in outgoing EF."""
    t = stress(r, **kw)
    f = f0(r,kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
    return dict(Trr=4*t['Tvv']/f**2,
                Tur=2*(t['Tvv']+t['Tuv'])/f,
                Tuu_ef=t['Tuu']+2*t['Tuv']+t['Tvv'])


def tortoise(r, r_reference=400., **kw):
    """Zero rstar at reference r0; classical tortoise in a static patch."""
    return quad(lambda s:1/f0(s,kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.)),
                r_reference,r,epsabs=1e-11,epsrel=1e-11)[0]


def outgoing_radial_solution(r,r_reference=400.,**kw):
    """sigma1^-(r), F1^-(u0,r) for sigma1(r0)=F1(u0,r0)=0."""
    G=kw.get('G',1.)
    sig=quad(lambda s:G*outgoing_stress(s,**kw)['Trr']/s,
             r_reference,r,epsabs=1e-12,epsrel=1e-11)[0]
    F1=quad(lambda s:-2*G*outgoing_stress(s,**kw)['Tur'],
            r_reference,r,epsabs=1e-12,epsrel=1e-11)[0]/r
    return sig,F1


def outgoing_radial_solution_ode(r,r_reference=400.,**kw):
    G=kw.get('G',1.)
    def rhs(s,y):
        T=outgoing_stress(s,**kw)
        return [G*T['Trr']/s,(-2*G*T['Tur']-y[1])/s]
    sol=solve_ivp(rhs,(r_reference,r),(0.,0.),rtol=2e-12,atol=1e-13)
    if not sol.success: raise RuntimeError(sol.message)
    return tuple(sol.y[:,-1])


def outgoing_F1(u,r,u_reference=0.,r_reference=400.,**kw):
    G=kw.get('G',1.)
    J=flux(kw.get('mass',100.),kw.get('L',3e-6),G,kw.get('N',12))
    return outgoing_radial_solution(r,r_reference,**kw)[1]+2*G*J*(u-u_reference)/r


def coordinate_match(v,r,v_reference=0.,r_reference=400.,**kw):
    """Same F=(grad r)^2 in two patches using u=v-2rstar +O(epsilon).

    F1 comparison at fixed spacetime event. Homogeneous mass term uses same
    normalization rF1(v0,r0)=rF1(u0,r0)=0.
    """
    u=v-2*tortoise(r,r_reference,**kw)
    ing=radial_solution(r,r_reference,**kw)[1]
    J=flux(kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.),kw.get('N',12))
    ing+=2*kw.get('G',1.)*J*(v-v_reference)/r
    out=outgoing_F1(u,r,v_reference,r_reference,**kw)
    return ing,out,u


def outgoing_field_residuals(r,**kw):
    """Check outgoing EF source identities in the convention E=2G*T."""
    f=f0(r,kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
    G=kw.get('G',1.)
    T=outgoing_stress(r,**kw)
    J=flux(kw.get('mass',100.),kw.get('L',3e-6),G,kw.get('N',12))
    # E_ur=-P1; E_uu=-f*P1+r*F1_u
    Er=2*G*T['Tur']
    Eu=2*G*f*T['Tur']+2*G*J
    return Er,2*G*T['Tur'],Eu,2*G*T['Tuu_ef']


if __name__=='__main__':
    print('r, sigma_in, sigma_out, F1_in, F1_out, F1_difference')
    for r in [260.,300.,400.,500.,650.,800.]:
        ing,out,u=coordinate_match(0.,r)
        print(r,radial_solution(r)[0],outgoing_radial_solution(r)[0],ing,out,ing-out)


def coordinate_correction(v,r,r_reference=400.,**kw):
    """First-order null-coordinate correction: u=v-2*rstar+chi+O(NG)^2.

    Homogeneous chart reparametrization fixed by chi(v, r_reference)=0.
    `v` is measured from reference time. Valid only in the overlap and for
    perturbatively small stress and accumulated secular mass loss.
    """
    def integrand(s):
        f=f0(s,kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
        sigma_minus=outgoing_radial_solution(s,r_reference,**kw)[0]
        Fplus=radial_solution(s,r_reference,**kw)[1]
        return 2*sigma_minus/f+2*Fplus/f**2
    c=quad(integrand,r_reference,r,epsabs=1e-12,epsrel=1e-11)[0]
    sigma_plus=radial_solution(r,r_reference,**kw)[0]
    sigma_minus=outgoing_radial_solution(r,r_reference,**kw)[0]
    return v*(sigma_plus-sigma_minus)+c


def integrability_identity(r,**kw):
    """Source-level check: sigma_+' - sigma_-'=4G J/(r*f^2)."""
    G=kw.get('G',1.)
    f=f0(r,kw.get('mass',100.),kw.get('L',3e-6),G)
    t=stress(r,**kw)
    delta_prime=G*(4*t['Tuu']/f**2-outgoing_stress(r,**kw)['Trr'])/r
    J=flux(kw.get('mass',100.),kw.get('L',3e-6),G,kw.get('N',12))
    return delta_prime,4*G*J/(r*f**2)
