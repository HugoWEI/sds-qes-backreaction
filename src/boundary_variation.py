"""Classical spherical SdS boundary variational audit.

S_bulk=1/(4G) integral sqrt(-g)[X R+(grad X)^2/(2X)+2(1-Lambda X)]
S_GHY=1/(2G) integral_boundary sqrt(-h) X K (outward normal).

In Gaussian normal ds²=dell²-A(t,ell)²dt², the GHY cancellation
leaves L_eff=(2 X_ell A_ell + A X_ell²/(2X)
            -X_t²/(2X A)+2A(1-Lambda X))/(4G).
A finite reference term S_ref=-1/(2G) integral dt A W(X) leads
to the on-shell radial boundary momenta below. Quasilocal E=-Pi_A
is a chosen classical reference-subtracted Brown-York energy. A fully
semiclassical energy must also include Polyakov boundary/auxiliary terms.
Not a gravity-to-bath junction and not a QES calculation.
"""
import numpy as np
from linearized_geometry import f0


def reference_W(X,L=3e-6):
    """2D boundary counterterm W=2r sqrt(f_dS), r=sqrt(X)."""
    if X <= 0 or L*X>=3:
        raise ValueError('reference worldtube must lie inside dS static horizon')
    return 2*np.sqrt(X)*np.sqrt(1-L*X/3)


def reference_W_prime(X,L=3e-6):
    y=1-L*X/3
    if X<=0 or y<=0:raise ValueError('invalid reference dilaton')
    return np.sqrt(y)/np.sqrt(X)-L*np.sqrt(X)/(3*np.sqrt(y))


def radial_boundary_momenta(A,X,X_n,A_n,G=1.,L=3e-6):
    """Radial canonical Pi_A and Pi_X after GHY and reference term.

    Radial Gaussian normal coordinate increases along the OUTWARD normal.
    Pi_A=(X_n-W)/(2G),
    Pi_X=A_n/(2G)+A X_n/(4G X)-A W_X/(2G).
    """
    if min(A,X,G)<=0:raise ValueError('positive A,X,G required')
    W=reference_W(X,L)
    wx=reference_W_prime(X,L)
    return ((X_n-W)/(2*G),
            A_n/(2*G)+A*X_n/(4*G*X)-A*wx/(2*G))


def classical_sds_boundary(b,mass=100.,L=3e-6,G=1.):
    f=f0(b,mass,L,G)
    if f<=0:raise ValueError('worldtube must be in SdS static patch')
    A=np.sqrt(f)
    fp=2*G*mass/b**2-2*L*b/3
    X=b*b;X_n=2*b*A;A_n=fp/2
    pA,pX=radial_boundary_momenta(A,X,X_n,A_n,G,L)
    return dict(lapse=A,Pi_A=pA,Pi_X=pX,quasilocal_energy=-pA,
                X_n=X_n,A_n=A_n,reference_W=reference_W(X,L))


def brown_york_reference_mass(b,mass=100.,L=3e-6,G=1.):
    """Reference-subtracted classical BY E; M=0 reference, same Lambda."""
    fref=1-L*b*b/3
    f=f0(b,mass,L,G)
    if min(f,fref)<=0:raise ValueError('outside common static patch')
    return b/G*(np.sqrt(fref)-np.sqrt(f))


def dEby_dM(b,mass=100.,L=3e-6,G=1.):
    """Exactly 1/sqrt(f), not including quantum Polyakov charges."""
    f=f0(b,mass,L,G)
    if f<=0:raise ValueError('outside static patch')
    return 1/np.sqrt(f)


def boundary_data_variations(A,X,X_n,A_n,G=1.,L=3e-6):
    """On-shell delta S=integral dt(Pi_A delta A+Pi_X delta X).

    With Dirichlet A and X both fixed, there is NO natural junction
    equation: an actual interface requires letting boundary data vary
    and specifying the nongravitating bath/controller contributions.
    """
    return radial_boundary_momenta(A,X,X_n,A_n,G,L)


if __name__=='__main__':
    for r in (300.,500.,750.):
        print(r,classical_sds_boundary(r),dEby_dM(r))
