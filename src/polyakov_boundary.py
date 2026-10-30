"""Boundary canonical data of the localized Polyakov action in 2D spherical SdS.

F(X,chi)=X/(4G)-N*chi/(48*pi);  S = ∫sqrt(-g)[F R + (∂X)^2/(8GX)
 +(1-Lambda X)/(2G) -N(∂chi)^2/(96*pi)] + 2∮sqrt(|h|) F K
 - (1/(2G))∮sqrt(|h|) W(X).

The auxiliary field is a LOCAL representation of the Polyakov anomaly: its
boundary condition/Green function is physical state input, not a new free
propagating particle to be given arbitrary independent data.

A Dirichlet-chi finite-radius boundary is mathematically well posed but may
not be a physical transparent interface. Natural free variation of chi requires
Pi_chi=0; the stipulated Unruh-de Sitter (UdS) future-horizon-regular state
instead has Pi_chi=-N*(kb-kc)/(48*pi), nonzero at nondegenerate SdS.

All charges refer to one chosen classical de Sitter subtraction W(X) and
Gaussian outward normal toward increasing r. No Polyakov finite state-dependent
counterterm, bath contribution, or self-consistent geometry is included.
"""
import numpy as np
from baseline import horizons, kappa, flux
from linearized_geometry import f0, f0_prime, f0_double_prime
from boundary_variation import reference_W, reference_W_prime, classical_sds_boundary


def coefficients(G=1., N=12):
    if G <= 0 or N <= 0:
        raise ValueError('positive G and central charge required')
    return dict(F_X=1/(4*G), F_chi=-N/(48*np.pi),kin_chi=-N/(96*np.pi))


def state_chi(t, r, mass=100., L=3e-6, G=1., N=12, ref_r=400.):
    """UdS choice chi=-log(f)-kb*u-kc*v, u=t-r*, v=t+r*.

    r* is normalized to 0 at ref_r; this shifts chi by a constant only.
    Valid strictly in the static patch. Signs make chi finite at each future
    horizon in its own EF chart, not at either past horizon.
    """
    from two_patch import tortoise
    f=f0(r,mass,L,G)
    if f<=0: raise ValueError('static patch only')
    kb,kc=kappa(mass,L,G)
    star=tortoise(r,ref_r,mass=mass,L=L,G=G,N=N)
    return -np.log(f)-(kb+kc)*t+(kb-kc)*star


def state_derivatives(r,mass=100.,L=3e-6,G=1.,N=12):
    f=f0(r,mass,L,G)
    if f<=0: raise ValueError('static patch only')
    kb,kc=kappa(mass,L,G)
    fp=f0_prime(r,mass,L,G)
    return dict(chi_t=-(kb+kc),chi_r=(kb-kc-fp)/f,
                chi_normal=(kb-kc-fp)/np.sqrt(f),
                future_bh_ingoing_r_slope=(2*kb-fp)/f,
                future_cosmo_outgoing_r_slope=(-2*kc-fp)/f)


def local_stress(r,mass=100.,L=3e-6,G=1.,N=12):
    """Polyakov chi reproduces chiral Tuu,Tvv,Tuv of the prescribed state."""
    f=f0(r,mass,L,G)
    if f<=0:raise ValueError('static patch only')
    fp=f0_prime(r,mass,L,G);fpp=f0_double_prime(r,mass,L,G)
    kb,kc=kappa(mass,L,G)
    tuu=N/(192*np.pi)*(2*f*fpp-fp*fp+4*kb*kb)
    tvv=N/(192*np.pi)*(2*f*fpp-fp*fp+4*kc*kc)
    tuv=N*f*fpp/(96*np.pi)
    scalar_curvature=-fpp
    trace=-4*tuv/f
    return dict(Tuu=tuu,Tvv=tvv,Tuv=tuv,trace=trace,
                anomaly=N*scalar_curvature/(24*np.pi),
                net_chiral_flux=tuu-tvv)


def canonical_momenta(A,X,X_n,A_n,chi_n,N=12,G=1.,L=3e-6):
    """Pi_(A,X,chi) of GHY-improved classical+Polyakov local action.

    Boundary finite reference W(X) does not depend on chi. A fixed dilaton
    and induced metric + fixed chi are Dirichlet; a free chi needs Pi_chi=0
    or a compensating bath/boundary functional.
    """
    if min(A,X,G,N)<=0:raise ValueError('positive A,X,G,N required')
    W=reference_W(X,L);Wx=reference_W_prime(X,L)
    return dict(Pi_A=(X_n-W)/(2*G)-N*chi_n/(24*np.pi),
                Pi_X=A_n/(2*G)+A*X_n/(4*G*X)-A*Wx/(2*G),
                Pi_chi=-N*(2*A_n+A*chi_n)/(48*np.pi))


def sds_boundary_data(b=500.,mass=100.,L=3e-6,G=1.,N=12):
    f=f0(b,mass,L,G)
    if f<=0:raise ValueError('worldtube outside static patch')
    A=np.sqrt(f);X=b*b
    der=state_derivatives(b,mass,L,G,N)
    p=canonical_momenta(A,X,2*b*A,f0_prime(b,mass,L,G)/2,
                        der['chi_normal'],N,G,L)
    kb,kc=kappa(mass,L,G)
    classical=classical_sds_boundary(b,mass,L,G)
    return dict(**p,**der,lapse=A,central_charge=N,
                classical_energy=classical['quasilocal_energy'],
                bare_polyakov_energy_correction=N*der['chi_normal']/(24*np.pi),
                bare_total_energy=-p['Pi_A'],
                natural_chi_boundary_residual=p['Pi_chi'],
                required_coordinate_chi_source=-p['Pi_chi'],
                required_proper_chi_source=-p['Pi_chi']/A,
                state_killing_flux=flux(mass,L,G,N),
                chi_boundary_time_derivative=-(kb+kc),
                chi_ward_term=p['Pi_chi']*(-kb-kc),
                kappa_b=kb,kappa_c=kc)


def horizon_regular_derivatives(which='bh',mass=100.,L=3e-6,G=1.,N=12):
    """Analytic EF derivative limits for the future-horizon-regular chi.

    BH ingoing: (∂r chi)_v=(2kb-f')/f -> -f''(rb)/f'(rb).
    Cosmo outgoing: (∂r chi)_u=(-2kc-f')/f -> -f''(rc)/f'(rc).
    """
    rb,rc=horizons(mass,L,G)
    r=rb if which=='bh' else rc if which=='cosmo' else None
    if r is None:raise ValueError('bh or cosmo')
    return -f0_double_prime(r,mass,L,G)/f0_prime(r,mass,L,G)


if __name__=='__main__':
    for r in (250.,350.,500.,700.,850.):
        d=sds_boundary_data(r)
        print('r',r,'Pi_chi',d['Pi_chi'],'Pi_chi*dotchi',d['chi_ward_term'],
              'flux',d['state_killing_flux'],
              'natural chi source (coordinate)',d['required_coordinate_chi_source'])
